# coding=utf-8

"""Produce realisations of the reaction-diffusion-master-equation

This algorithm is known as: kinetic Monte Carlo, stochastic simulation algorithm, Gillespie algorithm.
A proper description of such an algorithm is given by

    R. Erban, J. Chapman, and P. Maini, "A practical guide to stochastic simulations of reaction-diffusion processes", pp. 24–29, Apr. 2007.
"""

import bisect
import copy
import logging
import math
import tqdm as tqdm

import numpy as np

from pynumtools.kmc.logutil import StyleAdapter
from pynumtools.pynumtools_binding.kmc import convert_to_timeseries

__license__ = "LGPL"
__authors__ = ["chrisfroe", "clonker"]

log = StyleAdapter(logging.getLogger(__name__))


class Event:
    """Event is performed _on_ the system and applies a delta.

    The delta can be across boxes or across species, which is why its dimensions can be either (n_boxes,) or (n_species,).
    This behavior is defined by the derived classes.
    """

    def __init__(self, stoichiometric_delta, cumulative_propensity, interruptive=False):
        self._stoichiometric_delta = stoichiometric_delta
        self.cumulative_propensity = cumulative_propensity
        self.interruptive = interruptive

    def perform(self, state):
        pass

    def __repr__(self):
        return str(self)


class ReactionEvent(Event):
    def __init__(self, box_idx, stoichiometric_delta, cumulative_propensity, interruptive=False):
        super(ReactionEvent, self).__init__(stoichiometric_delta, cumulative_propensity, interruptive)
        self._box_idx = box_idx

    def __str__(self):
        string = "ReactionEvent\n"
        string += "box_idx " + str(self._box_idx) + "\n"
        string += "stoichiometric_delta (within one box, across species) " + str(self._stoichiometric_delta) + "\n"
        string += "cumulative_propensity " + str(self.cumulative_propensity)
        return string

    def perform(self, state):
        state[self._box_idx] += self._stoichiometric_delta


class DiffusionEvent(Event):
    def __init__(self, species_idx, stoichiometric_delta, cumulative_propensity, interruptive=False):
        super(DiffusionEvent, self).__init__(stoichiometric_delta, cumulative_propensity, interruptive)
        self._species_idx = species_idx

    def __str__(self):
        string = "DiffusionEvent\n"
        string += "species_idx " + str(self._species_idx) + "\n"
        string += "stoichiometric_delta (within one species, across boxes) " + str(self._stoichiometric_delta) + "\n"
        string += "cumulative_propensity " + str(self.cumulative_propensity)
        return string

    def perform(self, state):
        state[:, self._species_idx] += self._stoichiometric_delta


class Reaction:
    """Reactions behave like structs, that carry the rates and stoichiometric information for all possible reactions"""

    def __init__(self, rate, n_species, n_boxes, species_names):
        assert len(rate) == n_boxes
        self.rate = rate
        self.stoichiometric_delta = np.zeros(n_species, dtype=int)
        self.species_names = species_names

    def __repr__(self):
        return str(self)


class Conversion(Reaction):
    def __init__(self, species_from, species_to, rate, n_species, n_boxes, species_names):
        super(Conversion, self).__init__(rate, n_species, n_boxes, species_names)
        self.species_from = species_from
        self.species_to = species_to
        self.stoichiometric_delta[self.species_from] -= 1
        self.stoichiometric_delta[self.species_to] += 1

    def propensity(self, box_state, box_idx):
        return self.rate[box_idx] * box_state[self.species_from]

    def __str__(self):
        return "{} -> {}".format(self.species_names[self.species_from], self.species_names[self.species_to])

    def __repr__(self):
        string = "kmc.Conversion"
        string += "(from " + str(self.species_from) + ", to " + str(self.species_to) + ", rate " + str(self.rate) + ")"
        return string


class Fusion(Reaction):
    def __init__(self, species_from1, species_from2, species_to, rate, n_species, n_boxes, species_names):
        super(Fusion, self).__init__(rate, n_species, n_boxes, species_names)
        self.species_from1 = species_from1
        self.species_from2 = species_from2
        self.species_to = species_to
        self.stoichiometric_delta[self.species_from1] -= 1
        self.stoichiometric_delta[self.species_from2] -= 1
        if self.species_to is not None:
            self.stoichiometric_delta[self.species_to] += 1

    def propensity(self, box_state, box_idx):
        return self.rate[box_idx] * box_state[self.species_from1] * box_state[self.species_from2]

    def __str__(self):
        return "{} + {} -> {}".format(self.species_names[self.species_from1], self.species_names[self.species_from2],
                                      self.species_names[self.species_to] if self.species_to is not None else "0")

    def __repr__(self):
        string = "kmc.Fusion"
        string += "(from1 " + str(self.species_from1) + ", from2 " + str(self.species_from2)
        string += ", to " + str(self.species_to) + ", rate " + str(self.rate) + ")"
        return string


class Fission(Reaction):
    def __init__(self, species_from, species_to1, species_to2, rate, n_species, n_boxes, species_names):
        super(Fission, self).__init__(rate, n_species, n_boxes, species_names)
        self.species_from = species_from
        self.species_to1 = species_to1
        self.species_to2 = species_to2
        self.stoichiometric_delta[self.species_from] -= 1
        self.stoichiometric_delta[self.species_to1] += 1
        self.stoichiometric_delta[self.species_to2] += 1

    def propensity(self, box_state, box_idx):
        return self.rate[box_idx] * box_state[self.species_from]

    def __str__(self):
        return "{} -> {} + {}".format(self.species_names[self.species_from], self.species_names[self.species_to1],
                                      self.species_names[self.species_to2])

    def __repr__(self):
        string = "kmc.Fission"
        string += "(from " + str(self.species_from) + ", to1 " + str(self.species_to1)
        string += ", to2 " + str(self.species_to2) + ", rate " + str(self.rate) + ")"
        return string


class Decay(Reaction):
    def __init__(self, species_from, rate, n_species, n_boxes, species_names):
        super(Decay, self).__init__(rate, n_species, n_boxes, species_names)
        self.species_from = species_from
        self.stoichiometric_delta[self.species_from] -= 1

    def propensity(self, box_state, box_idx):
        return self.rate[box_idx] * box_state[self.species_from]

    def __str__(self):
        return "{} -> Ø".format(self.species_names[self.species_from])

    def __repr__(self):
        string = "kmc.Decay"
        string += "(from " + str(self.species_from) + ", rate " + str(self.rate) + ")"
        return string


class Creation(Reaction):
    def __init__(self, species_to, rate, n_species, n_boxes, species_names):
        super(Creation, self).__init__(rate, n_species, n_boxes, species_names)
        self.species_to = species_to
        self.stoichiometric_delta[self.species_to] += 1

    def propensity(self, box_state, box_idx):
        return self.rate[box_idx]

    def __str__(self):
        return "Ø -> {}".format(self.species_names[self.species_to])

    def __repr__(self):
        string = "kmc.Creation"
        string += "(to " + str(self.species_to) + ", rate " + str(self.rate) + ")"
        return string


class Fission2To3(Reaction):

    def __init__(self, species_from, species_to, rate, n_species, n_boxes, species_names):
        super().__init__(rate, n_species, n_boxes, species_names)
        assert isinstance(species_from, (list, tuple)) and len(species_from) == 2
        assert isinstance(species_to, (list, tuple)) and len(species_to) == 3
        self.species_from = species_from
        self.species_to = species_to
        for s in self.species_from:
            self.stoichiometric_delta[s] -= 1
        for s in self.species_to:
            self.stoichiometric_delta[s] += 1

    def propensity(self, box_state, box_idx):
        return self.rate[box_idx] * box_state[self.species_from[0]] * box_state[self.species_from[1]]

    def __str__(self) -> str:
        return f"{self.species_names[self.species_from[0]]} + {self.species_names[self.species_from[1]]} -> " \
            f"{self.species_names[self.species_to[0]]} + {self.species_names[self.species_to[1]]} " \
            f"+ {self.species_names[self.species_to[2]]}"

    def __repr__(self):
        string = "kmc.Fission2To3"
        string += "(from " + str(self.species_from) + ", to " + str(self.species_to) + ", rate " + str(self.rate) + ")"
        return string


class DoubleConversion(Reaction):

    def __init__(self, species_from, species_to, rate, n_species, n_boxes, species_names):
        super().__init__(rate, n_species, n_boxes, species_names)
        assert isinstance(species_from, (list, tuple)) and len(species_from) == 2
        assert isinstance(species_to, (list, tuple)) and len(species_to) == 2
        self.species_from = species_from
        self.species_to = species_to
        for s in self.species_from:
            self.stoichiometric_delta[s] -= 1
        for s in self.species_to:
            self.stoichiometric_delta[s] += 1

    def propensity(self, box_state, box_idx):
        return self.rate[box_idx] * box_state[self.species_from[0]] * box_state[self.species_from[1]]

    def __str__(self) -> str:
        return f"{self.species_names[self.species_from[0]]} + {self.species_names[self.species_from[1]]} -> " \
            f"{self.species_names[self.species_to[0]]} + {self.species_names[self.species_to[1]]} "

    def __repr__(self):
        string = "kmc.DoubleConversion"
        string += "(from " + str(self.species_from) + ", to " + str(self.species_to) + ", rate " + str(self.rate) + ")"
        return string


class Fusion3To2(Reaction):
    def __init__(self, species_from, species_to, rate, n_species, n_boxes, species_names):
        super().__init__(rate, n_species, n_boxes, species_names)
        assert isinstance(species_from, (list, tuple)) and len(species_from) == 3
        assert isinstance(species_to, (list, tuple)) and len(species_to) == 2
        self.species_from = species_from
        self.species_to = species_to
        for s in self.species_from:
            self.stoichiometric_delta[s] -= 1
        for s in self.species_to:
            self.stoichiometric_delta[s] += 1

    def propensity(self, box_state, box_idx):
        return self.rate[box_idx] * box_state[self.species_from[0]] * box_state[self.species_from[1]] \
               * box_state[self.species_from[2]]

    def __str__(self) -> str:
        return f"{self.species_names[self.species_from[0]]} + {self.species_names[self.species_from[1]]} + " \
            f"{self.species_names[self.species_to[2]]} -> {self.species_names[self.species_to[0]]} " \
            f"+ {self.species_names[self.species_to[1]]}"

    def __repr__(self):
        string = "kmc.Fusion3To2"
        string += "(from " + str(self.species_from) + ", to " + str(self.species_to) + ", rate " + str(self.rate) + ")"
        return string


class TrajectoryConfig(object):
    def __init__(self, types, reactions):
        self.types = types
        self.reactions = reactions

    def __str__(self):
        string = "TrajectoryConfig from kinetic_monte_carlo module\n"
        string += "  - particle types " + str(self.types) + "\n"
        string += "  - reactions\n"
        reactions_strings = map(str, self.reactions)
        for rs in reactions_strings:
            string += "      " + rs + "\n"
        return string


class ReactionDiffusionSystem:
    def __init__(self, diffusivity, n_species, n_boxes, init_state, init_time=0., species_names=None):
        diffusivity = np.asarray(diffusivity)
        init_state = np.asarray(init_state)

        assert n_species > 0
        assert n_boxes > 0
        # diffusivity can be a list of sparse matrices or a rank 3 tensor
        assert len(diffusivity) == n_species, "len diffusivity: {} != {}".format(n_species, len(diffusivity))
        assert diffusivity[0].shape == (n_boxes, n_boxes,)
        assert init_state.shape == (n_boxes, n_species)
        self._n_species = n_species
        self._n_boxes = n_boxes
        self._diffusivity = diffusivity
        self._reactions = []
        self._n_reactions = 0
        self._state = np.copy(init_state)
        self._init_state = np.copy(init_state)
        self._time = init_time
        self._event_list = []
        self._time_list = [init_time]
        self._state_list = [init_state]
        self._is_finalized = False
        if species_names:
            assert len(species_names) == n_species
        else:
            species_names = [str(x) for x in range(n_species)]
        self._species_names = species_names
        self._names_to_ids = dict()
        for idx, name in enumerate(self._species_names):
            self._names_to_ids[name] = idx

    def get_trajectory_config(self):
        """Generate a TrajectoryConfig object for interoperability with readdy_learn

        :return: TrajectoryConfig object
        """

        traj_config = TrajectoryConfig(self._names_to_ids, self._reactions)
        return traj_config

    def get_counts_config(self, n_frames=None, timestep=None):
        assert n_frames is not None or timestep is not None
        counts, times = self.convert_events_to_time_series2(time_step=timestep)
        # flatten out spatial dimension and convert to float
        counts = np.sum(counts, axis=1, dtype=float)
        config = self.get_trajectory_config()
        return counts, times, config

    def __str__(self):
        string = "ReactionDiffusionSystem\n"
        string += "--- n_species " + str(self._n_species) + "\n"
        string += "--- n_boxes " + str(self._n_boxes) + "\n"
        string += "--- diffusivity\n" + str(self._diffusivity) + "\n"
        string += "--- reactions\n" + str(self._reactions) + "\n"
        string += "--- state\n" + str(self._state) + "\n"
        # string += "--- init_state\n" + str(self._init_state) + "\n"
        string += "--- time " + str(self._time) + "\n"
        # string += "--- event_list\n" + str(self._event_list) + "\n"
        # string += "--- time_list\n" + str(self._time_list) + "\n"
        # string += "--- state_list\n" + str(self._state_list) + "\n"
        return string

    def _assure_not_finalized(self):
        if not self._is_finalized:
            return
        else:
            raise RuntimeError("System has already been finalized")

    def _id_from_name(self, type_name):
        return self._names_to_ids[type_name] if type_name is not None else None

    def add_conversion(self, species_from, species_to, rate):
        self._assure_not_finalized()
        species_from = self._id_from_name(species_from)
        species_to = self._id_from_name(species_to)
        conversion = Conversion(species_from, species_to, rate, self._n_species, self._n_boxes, self._species_names)
        self._reactions.append(conversion)
        self._n_reactions = len(self._reactions)

    def add_fusion(self, species_from1, species_from2, species_to, rate):
        self._assure_not_finalized()
        species_from1 = self._id_from_name(species_from1)
        species_from2 = self._id_from_name(species_from2)
        species_to = self._id_from_name(species_to)
        fusion = Fusion(species_from1, species_from2, species_to, rate, self._n_species, self._n_boxes,
                        self._species_names)
        self._reactions.append(fusion)
        self._n_reactions = len(self._reactions)

    def add_fission(self, species_from, species_to1, species_to2, rate):
        self._assure_not_finalized()
        species_from = self._id_from_name(species_from)
        species_to1 = self._id_from_name(species_to1)
        species_to2 = self._id_from_name(species_to2)
        fission = Fission(species_from, species_to1, species_to2, rate, self._n_species, self._n_boxes,
                          self._species_names)
        self._reactions.append(fission)
        self._n_reactions = len(self._reactions)

    def add_decay(self, species_from, rate):
        self._assure_not_finalized()
        species_from = self._id_from_name(species_from)
        decay = Decay(species_from, rate, self._n_species, self._n_boxes, self._species_names)
        self._reactions.append(decay)
        self._n_reactions = len(self._reactions)

    def add_creation(self, species_to, rate):
        self._assure_not_finalized()
        species_to = self._id_from_name(species_to)
        creation = Creation(species_to, rate, self._n_species, self._n_boxes, self._species_names)
        self._reactions.append(creation)
        self._n_reactions = len(self._reactions)

    def add_fission_2_to_3(self, species_from, species_to, rate):
        self._assure_not_finalized()
        species_from = [self._id_from_name(n) for n in species_from]
        species_to = [self._id_from_name(n) for n in species_to]
        reaction = Fission2To3(species_from, species_to, rate, self.n_species, self.n_boxes, self._species_names)
        self._reactions.append(reaction)
        self._n_reactions = len(self._reactions)

    def add_fusion_3_to_2(self, species_from, species_to, rate):
        self._assure_not_finalized()
        species_from = [self._id_from_name(n) for n in species_from]
        species_to = [self._id_from_name(n) for n in species_to]
        reaction = Fusion3To2(species_from, species_to, rate, self.n_species, self.n_boxes, self._species_names)
        self._reactions.append(reaction)
        self._n_reactions = len(self._reactions)

    def add_double_conversion(self, species_from, species_to, rate):
        self._assure_not_finalized()
        species_from = [self._id_from_name(n) for n in species_from]
        species_to = [self._id_from_name(n) for n in species_to]
        reaction = DoubleConversion(species_from, species_to, rate, self.n_species, self.n_boxes, self._species_names)
        self._reactions.append(reaction)
        self._n_reactions = len(self._reactions)

    @property
    def n_species(self):
        return self._n_species

    @property
    def n_boxes(self):
        return self._n_boxes

    @property
    def n_reactions(self):
        return self._n_reactions

    @property
    def reactions(self):
        return copy.deepcopy(self._reactions)

    @property
    def diffusivity(self):
        return copy.deepcopy(self._diffusivity)

    def simulate(self, n_steps=None, target_time=None, interruptive_reaction=lambda i, r: False,
                 interruptive_diff=lambda s, i, j: False):
        """
        Simulate up to n_steps events.

        If an interruptive event is performed, the simulation will stop beforehand
        and the current time will be returned. If successfully doing all steps, return None.
        """
        assert (n_steps is not None) ^ (target_time is not None)
        np.random.seed()

        log.info("Simulate for {} steps", n_steps)
        self._n_reactions = len(self._reactions)
        self._is_finalized = True

        t = 0

        interrupt = False
        while not interrupt:
            log.debug("Step {}, system {}", t, self)
            possible_events = []
            cumulative = 0.
            # gather reaction events
            for i in range(self._n_boxes):
                for r in range(self._n_reactions):
                    propensity = self._reactions[r].propensity(self._state[i], i)
                    if propensity > 0.:
                        cumulative += propensity
                        delta = self._reactions[r].stoichiometric_delta
                        shall_interrupt = interruptive_reaction(i, r)
                        possible_events.append(ReactionEvent(i, delta, cumulative, shall_interrupt))
            # gather diffusion events
            for s in range(self._n_species):
                for i in range(self._n_boxes):
                    for j in range(self._n_boxes):
                        if i != j:
                            propensity = self._diffusivity[s][i, j] * self._state[i, s]
                            if propensity > 0.:
                                cumulative += propensity
                                delta = np.zeros(self._n_boxes, dtype=int)
                                delta[i] = -1
                                delta[j] = +1
                                shall_interrupt = interruptive_diff(s, i, j)
                                possible_events.append(DiffusionEvent(s, delta, cumulative, shall_interrupt))

            if cumulative == 0.:
                log.info("No events possible / system is frustrated, at step {}", t)
                return

            # draw time and cumulative value
            event_time = (1. / cumulative) * np.log(1. / np.random.random())
            rnd = np.random.random() * cumulative

            # find event corresponding to rnd, that shall be performed
            cumulative_list = [x.cumulative_propensity for x in possible_events]
            event_idx = bisect.bisect_right(cumulative_list, rnd)
            event = possible_events[event_idx]
            log.debug("Cumulative list {}", cumulative_list)
            log.debug("Random number {}", rnd)
            log.debug("Event index for random number {}", event_idx)
            log.debug("Performing event {}", event)

            # save event and time to sequence
            self._event_list.append(event)
            self._time_list.append(self._time + event_time)

            # update system and save state
            self._time += event_time
            event.perform(self._state)
            self._state_list.append(np.copy(self._state))
            if event.interruptive:
                log.info("system interrupted after event {}, state is now {}", event, self._state)
                return self._time, event

            t += 1

            if target_time is not None:
                interrupt = self._time >= target_time
            else:
                interrupt = t >= n_steps
        return None

    def rollback(self, target_time):
        """Roll the system back to the event where event-time is smaller than the target_time"""
        for i in range(len(self._event_list) - 1, -1, -1):
            if self._time_list[i] < target_time:
                log.info("Rolling system back to {}, this corresponds to index {}", target_time, i)
                self._event_list = self._event_list[:i]  # no +1 since event is always one entry shorter than times
                self._time_list = self._time_list[:i + 1]
                self._time_list[-1] = target_time
                self._state_list = self._state_list[:i + 1]
                self._state = np.copy(self._state_list[-1])
                return

    @property
    def sequence(self):
        return self._event_list, self._time_list, self._state_list

    @property
    def smallest_dt(self):
        return np.min(np.diff(self._time_list))

    def convert_events_to_time_series2(self, time_step):
        n_frames = int(math.floor((self._time_list[-1] - self._time_list[0]) / time_step))
        result = np.zeros((n_frames, self._n_boxes, self._n_species), dtype=np.uint32)
        times = np.linspace(0, n_frames * time_step, num=n_frames, endpoint=False, dtype=float)

        times_list = np.array(self._time_list, dtype=float)
        state_list = np.array(self._state_list, dtype=np.uint32)
        convert_to_timeseries(result, times, times_list, state_list)

        # state = 0
        # nstates = len(self._state_list)
        # for ix in range(len(times)):
        #     t = times[ix]
        #     if t <= self._time_list[state]:
        #         result[ix] = self._state_list[state]
        #     else:
        #         while state < nstates and t > self._time_list[state]:
        #             state += 1
        #         result[ix] = self._state_list[state]

        return result, times

    def convert_events_to_time_series(self, time_step=None, n_frames=None):
        """
        Convert to an equidistant time series. Therefore either a time_step or
        a number of frames must be given.

        For usage with readdy_learn (w/o diffusion as of 2017-06-29), do: counts = np.sum(result, axis=1)

        :param time_step: the temporal difference between each frame
        :param n_frames: the number of frames
        :return: a trajectory of shape (n_frames, n_boxes, n_species)
        """
        if not ((time_step is not None) ^ (n_frames is not None)):
            raise RuntimeError("Either time_step (x)or n_frames must be given")
        if len(self._time_list) < 2:
            raise RuntimeError("Sample some events first")
        if n_frames:
            time_step = (self._time_list[-1] - self._time_list[0]) / float(n_frames)
        else:
            n_frames = math.ceil((self._time_list[-1] - self._time_list[0]) / time_step)
        log.info("got time step {} and n_frames {}".format(time_step, n_frames))

        result = np.zeros((n_frames, self._n_boxes, self._n_species), dtype=int)
        times = np.linspace(self._time_list[0], self._time_list[-1], n_frames)
        current_t = self._time_list[0]
        last_passed_event_time_idx = 0
        result[0, :, :] = self._init_state
        for t in range(1, n_frames):
            # increment the time and find out how many events we have passed
            current_t += time_step
            n_passed_events = 0
            if current_t > self._time_list[last_passed_event_time_idx + 1]:
                n_passed_events = 1
                # if we have skipped one event we might have skipped another
                while current_t > self._time_list[last_passed_event_time_idx + n_passed_events + 1]:
                    n_passed_events += 1
            last_passed_event_time_idx += n_passed_events
            result[t] = self._state_list[last_passed_event_time_idx]
        return result, times
