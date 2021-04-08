#include "ipf.hpp"

#include <iterator>
#include <fstream>

#include "smc.hpp"
#include "rng.hpp"
#include "eta.hpp" // time remaining
#include "macros.hpp" // RIGHT_HERE


void iterated_filtering(const Model & model,
      const std::vector<Timeseries> & txss,
      const std::vector<bool> & track_coffin,
      const AlgoParameters algpar,
      int threads, unsigned long seed) {
  // check some conditions
  if ( txss.size() != track_coffin.size() ) {
    throw std::invalid_argument("txss and track_coffin must have the same size" + RIGHT_HERE);
  }
  // init some variables and the RNG
  Rng rng(seed);
  int J = algpar.J;
  int M = algpar.M;
  int D = algpar.D; // duplicate final iterations
  int R = txss.size(); // the number of panels/subjects
  double cooldown_factor = 0.0;
  if ( M > 0 ) {
    cooldown_factor = pow(0.01, 1.0/M); // FIXME: 0.01 should be a parameter!
  }

  // open file for output (TODO: pass filename as argument)
  std::ofstream file;
  std::string filename = "data/out/ipf_result-" + model.getName() + ".xml";
  file.open(filename.c_str());
  if ( !file.good() ) {
    std::stringstream message;
    message << "unable to open file '" << filename << "'" << RIGHT_HERE;
    throw std::runtime_error(message.str());
  }
  file << "<iterated_filtering >" << std::endl;
  // names of free parameters
  file << "<initial_params >" << std::endl;
  file << model.par << std::endl;
  file << "</initial_params>" << std::endl;

  // create a particle filter object
  ParticleFilter pf(J, model, seed, threads);

  // create eta object to keep track of time
  EtaEstimator eta(M+D);

  // start iterated filtering
  for ( int m = 0; m < M+D; ++m ) {
    file << "<iterated_filtering_step >" << std::endl;

    std::list<std::pair<double, int>> cond_log_liks;

    // make a random permutation of [0, 1, ..., R-1]
    std::vector<int> permut(R);
    rng.shuffle(permut.data(), 0, R);

    // don't do a random walk at the last interation(s)
    if ( m >= M ) {
      pf.coolDown(0.0);
    }

    for ( int r : permut ) { // loop over the panel in random order
      file << "<subject id='" << txss[r].identifier << "' "
           << "idx='" << r << "' "
           << ">" << std::endl;
      int N = txss[r].size();
      /* make sure that the correct parameters are used
       * and lock all individual-level parameters except thr r-th
       */
      pf.selectIndivParams(r);
      // reset the state (potentially samples s0 using r-parameters)
      pf.reset(rng);
      // store the genealogy for the final iteration
      if ( m == M+D-1 ) {
        pf.storeGenealogy(algpar.dt); // dt is the time step for the paths
      }
      // choose the track_coffin mode
      pf.setCoffinTracker(track_coffin[r]);
      // run through timeseries of observations
      for ( int n = 0; n < N; ++n ) {
        std::cout << "m: " << m << " r: " << r << " n: " << n << std::endl;
        auto cond_log_lik = pf.step(txss[r].obs[n]);
        cond_log_liks.push_back(cond_log_lik);
        // print some info...
        std::cout << "log like: " << cond_log_lik.first
                  << "\t J valid: " << cond_log_lik.second
                  << "\t J eff: " << pf.getEffectiveJ()
                  << "\t J sim: " << pf.getEffectiveSampleSize()
                  << "\t J cof: " << pf.getJcoffin()
                  << std::endl;
        // write particle filter statistics to file
        file << pf << std::endl;
      } // n = 0, ..., N-1
      auto & Gen = pf.getGenealogy();
      if ( !Gen.empty() ) {
        for ( int g = 0; g < algpar.G; ++g ) {
          StatePath path = Gen.sampleStatePath(rng);
          file << "<path>" << std::endl;
          for ( State & s : path ) {
            file << s << std::endl;
          }
          file << "</path>" << std::endl;
        }
      }
      file << "</subject>" << std::endl;
    } // r = 0, ..., R-1

    // compute and print total log-likelihood
    double log_lik = 0.0;
    bool finite_log_lik = true;
    for ( auto & ll_j : cond_log_liks ) {
      if ( ll_j.second > 0 ) {
        log_lik += ll_j.first;
      } else {
        finite_log_lik = false;
      }
    }

    // print progress to screen
    eta.update();
    std::cout << "total log likelihood: " << log_lik << " "
              << (finite_log_lik ? " (valid)" : " (INVALID)") << " "
              << "ETA: " << eta << std::endl;

    // write to XML file
    file << "<log_lik "
         << "val='" << log_lik << "' "
         << "finite='" << std::boolalpha << finite_log_lik << std::noboolalpha << "' "
         << "/>";
    file << std::endl;

    /* print parameter range and median
     * TODO: print after every unit to track auto-correlation
     */
    file << "<param_range >" << std::endl;
    auto [par1, par2] = pf.getParamRange(ALPHA_RANGE);
    file << par1 << std::endl << par2 << std::endl;
    file << "</param_range >" << std::endl;
    // median
    file << "<param_median >" << std::endl;
    auto par = pf.getParamMedian();
    file << par << std::endl;
    file << "</param_median >" << std::endl;

    // closing tag for this iteration
    file << "</iterated_filtering_step>" << std::endl;

    /* reset happens at the start of replicate
     * just cool down the particle filter:
     */
    pf.coolDown(cooldown_factor);
  } // m = 1, ..., M

  file << "</iterated_filtering>" << std::endl;
  file.close();
}




void iterated_filtering(const Model & model,
      const std::vector<Timeseries> & txss,
      const AlgoParameters algpar,
      int threads, unsigned long seed) {
  // set all track_coffin bools to false by default
  std::vector<bool> track_coffin(txss.size(), false);
  // and use the function defined above
  iterated_filtering(model, txss, track_coffin, algpar, threads, seed);
}
