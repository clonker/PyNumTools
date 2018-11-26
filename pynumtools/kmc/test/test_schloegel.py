import numpy as np
import pynumtools.kmc as kmc


if __name__ == '__main__':
    if False:
        import matplotlib.pyplot as plt
        system = kmc.ReactionDiffusionSystem([[[0.]]], 1, 1, [[0]], species_names=["A"])
        system.add_fission_2_to_3(["A", "A"], ["A", "A", "A"], [0.18])
        system.add_fusion_3_to_2( ["A", "A", "A"], ["A", "A"], [2.5e-4])
        system.add_decay("A", [5])
        system.add_creation("A", [200])

        system.simulate(n_steps=1000000)

        counts, times, state = system.sequence
        state = np.array(state).squeeze()
        print(state.shape)
        plt.plot(times, state)
        plt.show()
