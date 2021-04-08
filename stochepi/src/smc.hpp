#ifndef SMC_HPP_
#define SMC_HPP_

#include <vector>
#include <utility>

#include "rng.hpp"
#include "parallelism.hpp"
#include "state.hpp"
#include "data.hpp"
#include "model.hpp"
#include "particle.hpp"
#include "genealogy.hpp"

// range for stats of the state variables and parameters
constexpr double ALPHA_RANGE = 0.95;


class ParticleFilter : public Printable {
public:
  ParticleFilter(int J, const Model & model, unsigned long seed, int threads);
  virtual ~ParticleFilter();
  ParticleFilter(const ParticleFilter & )=delete;
  /** Reset the state to the initial values,
   * but keeps the particles' Parameter values.
   */
  void reset(Rng & rng);
  /** Make sure that the right individual-level parameters are used.
   */
  void selectIndivParams(int r);
  /** Save the Genealogy and make paths with time interval dt.
   * Clears any existing Genealogy.
   */
  void storeGenealogy(double dt);
  /** extract the Genealogy. Might be empty.
   */
  const Genealogy & getGenealogy() const;
  /** Cool down the auxiliary random walk by a factor q.
   */
  void coolDown(double q);
  /** step returns a pair cointaing the log-likelihood estimate
   * and the number of particles that remained valid.
   * step evolves the particles forward in time and re-samples them
   * according to their weights.
   * @param obs the next observations
   * @param J the desired swarm size. If J=0, (default),
   *    use the previous swarm size.
   * @returns The function returns the LOG likelihood
   *    and the number of valid particles
   * if this number is ZERO, then the returned log likelihood
   *    is invalid (because -inf)
   */
  std::pair<double, int> step(const Observation & obs, int J=0);
  // properties
  int getJ() const; // particles.size()
  int getEffectiveJ() const; // effective_swarm_size
  double getEffectiveSampleSize() const;
  int getJcoffin() const;
  /** pass true to treat the coffin state separately.
   * For instance, this avoids particle depletion in the case
   * where it is known that the coffin state becomes invalid
   * at some point.
   */
  void setCoffinTracker(bool track_coffin);
  void print(std::ostream & os) const override;
  void printAllParticles(std::ostream & os) const;
  /** Compute the quantiles of the current particle swarm.
   */
  std::vector<State> getQuantiles(const std::vector<double> & qs) const;
  // specialized functions for range and median (call getQuantiles)
  std::pair<State, State> getRange(double alpha) const;
  State getMedian() const;
  /** Compute the quantiles of parameters of the swarm
   */
  std::pair<Parameters, Parameters> getParamRange(double alpha) const;
  Parameters getParamMedian() const;
protected:
  WorkerPool wp;
  Model model;
  std::vector<Particle*> particles;
  int J_coffin; // special "extinct" particles
  bool track_coffin; // treat coffin particles differently
  Rng rng;
  double rel_temp; // for cooling down
  int effective_swarm_size; // number of unique particles
  double effective_sample_size; // inverse simpson index of the weights
  double cond_log_lik_hat; // log swarm-average conditional likelihood
  // if we don't want to print all particles, set print_modulus > 1
  int print_modulus;
  // flag for keeping track of the Particle Genealogy
  bool store_genealogy;
  // if store_genealogy is true, the Genealogy G will contain the filtered paths
  Genealogy G;
  // keep track of statistics: prediction range, median
  std::pair<State, State> pred_range;
  State pred_median;
  State filter_median;
};

// return the same number of indices as the number of given weights
std::map<int, int> systematicResampling(const std::map<int, double> & weights, Rng & rng);
// return a given number of weights J
std::map<int, int> systematicResampling(const std::map<int, double> & weights, int J, Rng & rng);

// return the same number of indices as the number of given weights (pass log weights)
std::map<int, int> systematicResamplingLogScale(const std::map<int, double> & logweights, Rng & rng);
// return a given number of weights J
std::map<int, int> systematicResamplingLogScale(const std::map<int, double> & logweights, int J, Rng & rng);



#endif
