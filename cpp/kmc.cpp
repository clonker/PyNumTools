/**
 * @file kmc.cpp
 * @brief Functions used in Kinetic Monte Carlo algorithm
 * @author chrisfroe
 * @date 28.11.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include "kmc.h"

namespace kmc {

void convert_to_timeseries(kmc_result_array &result, const kmc_times_array &times,
                                  const kmc_times_array &times_list, const kmc_state_array &state_list) {
    std::size_t state = 0;

    auto nstates = (std::size_t) state_list.shape()[0];
    auto nframes = (std::size_t) result.shape()[0];
    auto nboxes = (std::size_t) result.shape()[1];
    auto nspecies = (std::size_t) result.shape()[2];

    for(std::size_t ix = 0; ix < nframes; ++ix) {
        auto t = times.at(ix);
        if(t <= times_list.at(state)) {
            for (std::size_t s = 0; s < nspecies; ++s) {
                for (std::size_t b = 0; b < nboxes; ++b) {
                    result.mutable_at(ix, b, s) = state_list.at(state, b, s);
                }
            }
        } else {
            while(state < nstates && t > times_list.at(state)) {
                ++state;
            }
            for (std::size_t s = 0; s < nspecies; ++s) {
                for (std::size_t b = 0; b < nboxes; ++b) {
                    result.mutable_at(ix, b, s) = state_list.at(state, b, s);
                }
            }
        }
    }
}

}