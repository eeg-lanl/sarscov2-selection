#ifndef PARTICLE_HPP_
#define PARTICLE_HPP_

#include <vector>
#include <utility>
#include <mutex>

#include "rng.hpp"
#include "parallelism.hpp"
#include "state.hpp"
#include "engine.hpp"
#include "data.hpp"
#include "model.hpp"

typedef std::list<State> StatePath;

typedef long int Id;
constexpr Id INVALID_ID = -1;

class Particle : public Printable {
public:
  Particle(const State & s, const Parameters & par, const Engine & e,
    const ObservableFunction & obsfun, bool sde);
  /** We must be able to copy particles, but a Job object cannot be copied
   * due to mutex members and hence the default copy constructor
   * is deleted by default
   */
  Particle(const Particle & particle);
  Particle & operator=(const Particle & particle);
  virtual ~Particle() { /* empty */ }
  const State & getState() const;
  State & getState();
  const State & getObservables() const;
  const Parameters & getParameters() const;
  const double & getCachedLogPrior() const;
  const std::pair<double, bool> & getCachedLogLike() const;
  void cacheLogPrior(double );
  void cacheLogLike(const std::pair<double, bool> & );
  void setNextTime(double t_next); // update t_next
  void setRelTemp(double rel_temp); // update rel_temp
  /** set store_path to true and dt_path to dt.
   * Throws exception in dt <= 0.
   */
  void storePath(double dt);
  /// Clear the stored path and sets store_path to false
  void clearPath();
  /** extract the path. TODO: method that requires less copying
   */
  StatePath getPath() const;
  void mutatePar(Rng & rng);
  void selectIndivParams(int r);
  void advanceState(Rng & rng);
  void resetState(const State & s); // also resets t_next to s.t()
  void getNewId(); // make id_parent equal to id and use next_id to get a new ID
  void print(std::ostream & os) const override;
  // execute is used by smc to mutate and advance...
  void execute(Rng & rng);
protected:
  /** auxiliary method used by copy constructor
   * and copy assignment constructor.
   * Copies only members with a default constructor.
   */
  void copy(const Particle & particle);
  // members
  State s;
  Parameters par;
  Engine e;
  ObservableFunction obsfun;
  State obs;
  double t_next; // the destination time of execute
  double rel_temp; // the relative temperature used for
  Id id; // this particle's id
  Id id_parent; // the id of the parent
  std::pair<double, bool> loglike_cache;
  double logprior_cache;
  /** For visualization, the Particle can remember the path leading to it.
   * store_path determines if this happens.
   */
  bool store_path;
  /** The sampling interval of the optional path is given by dt_path */
  double dt_path;
  /** The actual States are stored in a list */
  StatePath path;
  /** what kind of advancer is used? */
  bool sde;
  /** id counter shared between classes, guarded by a mutex */
  static std::mutex next_id_mutex;
  volatile static Id next_id;
};


#endif
