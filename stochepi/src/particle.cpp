#include "particle.hpp"


std::mutex Particle::next_id_mutex;
volatile Id Particle::next_id = 0; // the static Particle counter

// METHODS FOR Particle

Particle::Particle(const State & s, const Parameters & par,
    const Engine & e, const ObservableFunction & obsfun,
    bool sde) : s(s), par(par), e(e), obsfun(obsfun), sde(sde) {
  // compute observables using obsfun, state and parameters
  obs = obsfun(this->s, this->par);
  t_next = this->s.t();
  rel_temp = 1.0;
  id_parent = INVALID_ID;
  /* set the cached loglike and logprior to a VALID 0.0,
   * such that we can use it in the first step of iterated filtering
   */
  logprior_cache = 0.0;
  loglike_cache = std::make_pair(0.0, true);
  // by default: do not store any intermediate states
  store_path = false;
  dt_path = 0.0;
  // get the next Id, but first lock the corresponding mutex
  std::lock_guard<std::mutex> lock(next_id_mutex);
  id = next_id++;
}

Particle::Particle(const Particle & particle) : s(particle.s),
    par(particle.par), e(particle.e),
    obsfun(particle.obsfun), obs(particle.obs) {
  copy(particle);
}

Particle & Particle::operator=(const Particle & particle) {
  if ( this != &particle ) {
    s = particle.s;
    par = particle.par;
    e = particle.e;
    obsfun = particle.obsfun;
    obs = particle.obs;
    copy(particle);
  } // else do nothing
  return *this;
}

void Particle::copy(const Particle & particle) {
  // only copy members with default copy constructors
  t_next = particle.t_next;
  rel_temp = particle.rel_temp;
  id_parent = particle.id;
  logprior_cache = particle.logprior_cache;
  loglike_cache = particle.loglike_cache;
  store_path = particle.store_path;
  dt_path = particle.dt_path;
  path = particle.path;
  sde = particle.sde;
  // get the next Id, but first lock the corresponding mutex
  std::lock_guard<std::mutex> lock(next_id_mutex);
  id = next_id++;
}

void Particle::mutatePar(Rng & rng) {
  par.mutate(rng, rel_temp, s.t());
}

void Particle::selectIndivParams(int r) {
  // make sure that par returns the correct value
  par.select(r);
  // unlock all inidividual parameters except the focal one
  par.lockAllBut(r);
}

void Particle::advanceState(Rng & rng) {
  if ( s.t() < t_next ) {
    if ( store_path ) {
      // remove old path
      path.clear();
      // save the current state
      path.push_back(concat(s, obs));
      while ( s.t() + dt_path < t_next ) {
        double t_inter = s.t() + dt_path;
        s = e.evolve_hybrid(s, t_inter, par, rng, sde);
        State obs_inter = obsfun(s, par);
        path.push_back(concat(s, obs_inter));
      }
    } // in any case: take the step to t_next
    s = e.evolve_hybrid(s, t_next, par, rng, sde);
    // compute observables using advanced state
    obs = obsfun(s, par);
  } // else do nothing
}

void Particle::resetState(const State & s) {
  this->s = s;
  t_next = s.t();
  e.reset(s, par);
  // reset observables
  obs = obsfun(this->s, par);
}


void Particle::execute(Rng & rng) {
  // NB: FIRST MUTATE PARAMETERS
  mutatePar(rng);
  // THEN EVOLVE THE STATE
  advanceState(rng);
}

const State & Particle::getState() const {
  return s;
}

State & Particle::getState() {
  return s;
}

const State & Particle::getObservables() const {
  return obs;
}

StatePath Particle::getPath() const {
  if ( path.empty() ) {
    // return a point-path
    std::list<State> point;
    point.push_back(concat(s, obs)); // TODO concatenate observables
    return point;
  } else {
    return path;
  }
}

const Parameters & Particle::getParameters() const {
  return par;
}

const double & Particle::getCachedLogPrior() const {
  return logprior_cache;
}
const std::pair<double, bool> & Particle::getCachedLogLike() const {
  return loglike_cache;
}
void Particle::cacheLogPrior(const double logprior) {
  logprior_cache = logprior;
}
void Particle::cacheLogLike(const std::pair<double, bool> & loglike) {
  loglike_cache = loglike;
}



void Particle::setNextTime(double t_next) { // update t_next
  this->t_next = t_next;
}

void Particle::setRelTemp(double rel_temp) { // update the relative temperature
  this->rel_temp = rel_temp;
}

void Particle::storePath(double dt) {
  store_path = true;
  if ( dt <= 0 ) {
    throw std::invalid_argument("path sampling interval must be positive");
  }
  dt_path = dt;
}

void Particle::clearPath() {
  store_path = false;
  dt_path = 0.0;
  path.clear();
}


void Particle::getNewId() {
  id_parent = id;
  // get the next Id, but first lock the corresponding mutex
  std::lock_guard<std::mutex> lock(next_id_mutex);
  id = next_id++;
}

void Particle::print(std::ostream & os) const {
  os << "<particle "
     << "id='" << id << "' "
     << "id_parent='" << id_parent << "' "
     << ">" << std::endl;
  // print state
  os << concat(s, obs) << std::endl;
  // print parameter values
  os << par << std::endl;
  // print path if stored
  if ( store_path ) {
    os << "<path "
       << "dt='" << dt_path << "' "
       << '>' << std::endl;
    for ( auto & state : path ) {
      os << state << std::endl;
    }
    os << "</path>" << std::endl;
  }
  // close xml element
  os << "</particle>";
}
