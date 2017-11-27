# coding=utf-8

import unittest
import numpy as np

import pynumtools.kmc as kmc

__license__ = "LGPL"
__authors__ = ["chrisfroe", "clonker"]

def example_system_conversions():
    """Produce realisations of system with only conversions"""
    n_species = 2
    n_boxes = 2
    diffusivity_0 = np.array([[0., 0.3], [0.4, 0.]])  # species 0
    diffusivity_1 = np.array([[0., 0.5], [0.9, 0.]])  # species 1
    diffusivity = np.array([diffusivity_0, diffusivity_1])
    init_state = np.array([[1, 1], [2, 2]], dtype=np.int)
    species_names = ["A", "B"]
    system = kmc.ReactionDiffusionSystem(diffusivity, n_species, n_boxes, init_state, species_names=species_names)
    system.add_conversion("A", "B", np.array([4., 4.]))
    system.add_conversion("B", "A", np.array([0.5, 0.5]))
    system.simulate(50)
    return system


def set_up_system(init_state):
    desired_rates = np.array([2e-2, 2e-2, 1e-2, 1e-4, 1e-4, 0,0,0,0,0,0,0,0,0,0, 0,0,0])

    sys = kmc.ReactionDiffusionSystem(n_species=4, n_boxes=1, diffusivity=[[[0.]], [[0.]], [[0.]], [[0.]]],
                                      init_state=init_state, species_names=["A", "B", "C", "D"])
    sys.add_conversion("A", "D", np.array([desired_rates[0]]))
    sys.add_conversion("D", "A", np.array([desired_rates[1]]))
    sys.add_conversion("D", "B", np.array([desired_rates[2]]))
    sys.add_fusion("A", "B", "C", np.array([desired_rates[3]]))
    sys.add_fission("C", "D", "B", np.array([desired_rates[4]]))

    return sys


class TestKineticMonteCarlo(unittest.TestCase):
    def test_raise_if_finalized(self):
        with np.testing.assert_raises(RuntimeError):
            n_species, n_boxes = 2, 2
            diffusivity_0 = np.array([[0., 0.3], [0.4, 0.]])  # species 0
            diffusivity_1 = np.array([[0., 0.5], [0.9, 0.]])  # species 1
            diffusivity = np.array([diffusivity_0, diffusivity_1])
            init_state = np.array([[1, 1], [2, 2]], dtype=np.int)
            system = kmc.ReactionDiffusionSystem(diffusivity, n_species, n_boxes, init_state)
            system.simulate(5)
            system.add_creation("1", 1.)

    def test_always_positive_number_of_particles(self):
        event_list, time_list, state_list = example_system_conversions().sequence
        state_array = np.asarray(state_list, dtype=np.int)
        all_positive = state_array >= 0
        self.assertTrue(np.all(all_positive))

    def test_conservation_of_particles(self):
        """In the system of only Conversion reactions, the total number of particles is conserved."""
        event_list, time_list, state_list = example_system_conversions().sequence
        n_particles = np.sum(state_list[0])
        all_correct = np.fromiter(map(lambda state: np.sum(state) == n_particles, state_list), dtype=np.bool)
        self.assertTrue(np.all(all_correct))

    def test_convert_to_time_series_args(self):
        with np.testing.assert_raises(Exception):
            system = example_system_conversions()
            system.convert_events_to_time_series([])

        with np.testing.assert_raises(Exception):
            system = example_system_conversions()
            system.convert_events_to_time_series([], time_step=0.1, n_frames=10)

    def test_conservation_of_particles_after_converting(self):
        system = example_system_conversions()
        time_series, times = system.convert_events_to_time_series(n_frames=500)
        n_particles = np.sum(time_series[0])
        all_correct = np.fromiter(map(lambda state: np.sum(state) == n_particles, time_series), dtype=np.bool)
        self.assertTrue(np.all(all_correct))

    def test_rollback(self):
        system = example_system_conversions()
        events, times, states = system.sequence
        target = times[-3] - 0.000001
        system.rollback(target)
        new_events, new_times, new_states = system.sequence
        self.assertEqual(len(events) - len(new_events), 3)
        self.assertEqual(len(times) - len(new_times), 3)
        self.assertEqual(len(states) - len(new_states), 3)
        np.testing.assert_equal(events[-4], new_events[-1])
        np.testing.assert_equal(target, new_times[-1])
        np.testing.assert_equal(states[-4], new_states[-1])

    def test_interrupt(self):
        system = example_system_conversions()
        init_state = np.copy(system._state)
        # interrupt after any reaction in box 0 or any diffusion from or to box 0
        # i.e. population in box 0 does change exactly once
        res = system.simulate(100, lambda i, r: i == 0, lambda s, i, j: i == 0 or j == 0)
        state = np.copy(system._state)
        self.assertTrue(res is not None)
        a_difference = abs(state[0][0] - init_state[0][0])
        b_difference = abs(state[0][1] - init_state[0][1])
        # A or B (or both ) have a delta of 1 in box 0
        self.assertTrue(a_difference == 1 or b_difference == 1)


if __name__ == '__main__':
    unittest.main()