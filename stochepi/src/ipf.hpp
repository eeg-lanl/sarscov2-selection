#ifndef IPF_HPP_
#define IPF_HPP_

#include "model.hpp"
#include "data.hpp"

/** A structure with some algorithmic parameters
 */
struct AlgoParameters {
      int J; // number of particles
      int M; // number of iterations
      int D; // number of duplicate runs after complete cooldown
      int G; // number of sampled trajectories
      double dt; // time step for the final PF run: store the genealogy!
};

void iterated_filtering(const Model & model,
      const std::vector<Timeseries> & txss,
      const std::vector<bool> & track_coffin,
      const AlgoParameters algpar,
      int threads, unsigned long seed);

/** set track_coffin to false */
void iterated_filtering(const Model & model,
      const std::vector<Timeseries> & txss,
      const AlgoParameters algpar,
      int threads, unsigned long seed);



#endif
