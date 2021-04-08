#include "engine.hpp"

#include <algorithm>
#include <numeric>
#include <cmath>
#include <functional> // for std::cref (pass reference of VectorField to odeint)
#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/generation.hpp> // make_controlled

#include "variable.hpp"
#include "sde.hpp" // SDE integrator

namespace odeint = boost::numeric::odeint;

typedef odeint::runge_kutta_cash_karp54<RealVec> real_vec_error_stepper_type;

constexpr bool VERBOSE = false; // TODO: better way to handle message printing

const double Engine::dt_max = 1.0;
const double Engine::dt_init = 0.001;
const double Engine::h_init = 0.01; // TODO: can be relaxed with higher-order solvers!

Engine::Engine(const State & s0, const std::vector<Transition*> transitions,
    const Parameters & par) : s0(s0), h(h_init), dt(dt_init), vf(s0, transitions, par) {
  // duplicate transitions
  replace_transitions(transitions); // copy transitions to this->transitions
}

Engine::Engine(const Engine & e) : s0(e.s0), h(e.h), dt(e.dt), vf(e.vf) {
  // copy constructor: duplicate transitions
  replace_transitions(e.transitions);
}

Engine & Engine::operator=(const Engine & e) { // copy assignment constructor
  if ( this != &e ) {
    s0 = e.s0;
    h = e.h;
    dt = e.dt;
    vf = e.vf;
    replace_transitions(e.transitions);
  } // else do nothing
  return *this;
}

void Engine::reset(const State & s0, const Parameters & par) {
  this->s0 = s0;
  // FIXME: reset algorithmic parameters??
  // FIXME: update VectorField?
}


State Engine::evolve_gillespie(State s, double tmax, const Parameters & p, Rng & rng) {
  while ( true ) {
    // choose the next event
    std::vector<double> rates(transitions.size(), 0.0);
    std::transform(transitions.begin(), transitions.end(), rates.begin(),
        [&] (Transition* x) -> double {return x->rate(s, p);});
    double lambda = std::accumulate(rates.begin(), rates.end(), 0.0);
    // check that lambda > 0 (otherwise jump to tmax)
    if ( lambda == 0 ) {
      s.t() = tmax;
      break;
    } // else, lambda > 0
    double time_increment = rng.exponential(lambda);
    // check that the next event takes place before tmax
    if ( std::nextafter(s.t() + time_increment, tmax) >= tmax ) {
      s.t() = tmax;
      break;
    } // else, apply the transition and incerment the time with dt

    Transition* trans = choose_random_transition(transitions, rates, lambda, rng);
    if ( trans == nullptr ) {
      throw std::logic_error("unable to randomly choose transition" + RIGHT_HERE);
    }
    s = trans->apply(s);
    s.t() += time_increment; // manually increase the time
    // TODO: make dt optional argument for Transition::apply?
  }
  return s;
}

State Engine::evolve_tauleap(State s, double tmax, const Parameters & p, Rng & rng) {
  while ( std::nextafter(s.t(), tmax) < tmax ) {
    // choose the next event
    std::vector<double> rates(transitions.size(), 0.0);
    std::transform(transitions.begin(), transitions.end(), rates.begin(),
        [&] (Transition* x) -> double {return x->rate(s, p);});
    double lambda = std::accumulate(rates.begin(), rates.end(), 0.0);
    // sample number of events
    int n = rng.poisson(lambda * dt);
    // do n transitions
    for ( int i = 0; i < n; ++i ) {
      Transition* trans = choose_random_transition(transitions, rates, lambda, rng);
      if ( trans == nullptr ) {
        throw std::logic_error("unable to randomly choose transition" + RIGHT_HERE);
      }
      s = trans->apply(s);
      // re-compute the rates..
      std::transform(transitions.begin(), transitions.end(), rates.begin(),
          [&] (Transition* x) -> double {return x->rate(s, p);});
      lambda = std::accumulate(rates.begin(), rates.end(), 0.0);
    }
    s.t() += dt; // manually increase the time
    // adjust the timestep
    if ( lambda > 0 ) {
      dt = std::min({dt_max, 1/lambda, tmax - s.t()}); // TODO: smoothen, use intermediate rates
    } else {
      dt = std::min(dt_max, tmax - s.t());
    }
  }
  return s;
}

State Engine::evolve_ode(State s, double tmax, const Parameters & p) {
  // make sure that all variables of the state are continuous
  for ( auto & var : s ) {
    var.make_continuous(); // make_continous makes the entire vector continuous
  }
  // initialize a vectorfield with the fully continous state
  vf.update(s, p);
  // integrate the system of ODEs
  if ( tmax > s.t() ) {
    std::tie(s, std::ignore) = integrate(tmax - s.t(), vf, s, h);
  }
  return s;
}


State Engine::evolve_sde(State s, double tmax, const Parameters & p, Rng & rng) {
  // make sure that all variables of the state are continuous
  for ( auto & var : s ) {
    var.make_continuous(); // make_continous makes the entire vector continuous
  }
  // initialize a vectorfield with the fully continous state
  vf.update(s, p);
  // integrate the system of SDEs
  if ( tmax > s.t() ) {
    std::tie(s, std::ignore) = sde_integrate(tmax - s.t(), vf, s, h, rng);
  }
  return s;
}

State Engine::evolve_hybrid(State s, double tmax,
    const Parameters & p, Rng & rng, bool sde) {
  // use THRESHOLD to define a backwards-compatible switcher
  VarTypeSwitcher switcher = [](const State & s, const Parameters & p) -> std::vector<bool> {
    std::vector<bool> ought_discrete(s.flat_size());
    for ( size_t i = 0; i < s.flat_size(); ++i ) {
      const Variable & var = s[i];
      if ( var.is_discrete )
        // keep discrete if below threshold + window/2
        ought_discrete[i] = (var.disc_value < THRESHOLD + WINDOW/2);
      else { // continuous
        // make discrete if below threshold - window/2
        ought_discrete[i] = (var.cont_value < THRESHOLD - WINDOW/2);
      }
    }
    return ought_discrete;
  };
  // now use the general method
  return evolve_hybrid(s, tmax, p, switcher, rng, sde);
}



State Engine::evolve_hybrid(State s, double tmax,
    const Parameters & p, const VarTypeSwitcher & switcher, Rng & rng, bool sde) {
  // make sure that the timestep "fits"
  while ( std::nextafter(s.t(), tmax) < tmax ) {
    /* test if we need to switch Variables from discrete
     * to continuous and vice versa. Do this first to
     * avoid problems with the initially defined state.
     */
    std::vector<bool> ought_discrete = switcher(s, p);
    for ( size_t i = 0; i < s.flat_size(); ++i ) {
      Variable & var = s[i];
      if ( var.is_discrete && !ought_discrete[i] ) {
        var.make_continuous();
        // print message for debugging etc.
        if ( VERBOSE ) {
          std::cout << "# switching from discrete to continuous" << std::endl;
        }
      }
      if ( !var.is_discrete && ought_discrete[i] ) {
        if ( VERBOSE && var.value() < 0 ) {
          std::cerr << "# WARNING: negative var detected: " << var << std::endl;
        }
        var.make_discrete();
        // be careful with the initial timestep...
        dt = std::min(dt, dt_init);
        // print message for debugging etc.
        if ( VERBOSE ) {
          std::cout << "# switching from continuous to discrete" << std::endl;
        }
      }
    }
    // re-assign indices whenever the signature of s has changed
    vf.update(s, p);
    // take a good time step based on instantaneous transition rate
    double lambda = vf.totalTransitionRate();
    if ( lambda > 0 ) {
      dt = std::min(dt_max, 1/lambda);
    } else {
      dt = dt_max;
    }
    double time_increment = std::min(dt, tmax - s.t());
    // integrate the deterministic variables and the loads
    std::map<Transition*, double> loads;
    if ( sde ) {
      // h is passed by ref and modified (TODO adaptive SDE integrator)
      std::tie(s, loads) = sde_integrate(time_increment, vf, s, h, rng);
    } else { // use ODE integrator
      // h is passed by ref and modified
      std::tie(s, loads) = integrate(time_increment, vf, s, h);
    }
    auto op = [] (double c, std::pair<Transition*, double> l) {return c + l.second;};
    double Lambda = std::accumulate(loads.begin(), loads.end(), 0.0, op);
    // sample number of events (0 if Lambda = 0)
    int n = (Lambda > 0 ? rng.poisson(Lambda) : 0); // Lambda = lambda * dt

    // do n transitions
    for ( int i = 0; i < n; ++i ) {
      if ( i > 0 ) {
        // re-compute the rates in order to apply another transition
        for ( auto & [trans, load] : loads ) {
          load = trans->rate(s, p) * dt; // approximate load by rate * dt
        }
        Lambda = std::accumulate(loads.begin(), loads.end(), 0.0, op);
        if ( Lambda == 0 ) {
          // print message for debugging etc.
          if ( VERBOSE ) {
            std::cerr << "# WARNING: sampling more than one event is not possible"
                      << RIGHT_HERE << std::endl;
          }
          break; // break for ( i = 0; i < n; ++i )
        }
      }
      // sample a random transition
      Transition* trans = choose_random_transition(loads, Lambda, rng);
      if ( trans == nullptr ) {
        throw std::logic_error("unable to choose random transition" + RIGHT_HERE);
      } // else..
      State s_prime = trans->apply(s);
      if ( s_prime.isNonNegative() ) {
        s = s_prime;
      } else {
        throw std::logic_error("about to do an illegal transition" + RIGHT_HERE);
      }
    } // for i = 0 ... n-1
  }
  return s;
}

/** Auxiliary integrate function is a wrapper around the Boost odeint methods.
 *
 * The State is "encoded" into a single vector, that can be handled by
 * odeint. After integration, the vector is decoded back into the State
 * object.
 */
std::pair<State, std::map<Transition*, double>>
    integrate(double dt, const VectorField & vf, const State & s, double & h) {
  // TODO: verify that s and vf are compatible
  auto controlled_stepper = odeint::make_controlled(1e-6, 1e-6, real_vec_error_stepper_type());
  // TODO: make the stepper member of engine
  RealVec y = vf.encode(s);

  double t_start = s.t();
  double t_end = s.t() + dt;
  double t_old = s.t(); // used to keep track of the timestep
  double h0 = h; // the initially proposed stepsize
  // observer is used to get h before final step to t_end
  auto observer = [&] (const RealVec & y, double t_new) {
    if ( t_old < t_new && t_new < t_end ) {
      h = t_new - t_old;
      t_old = t_new;
    }
  };
  size_t steps = odeint::integrate_adaptive(controlled_stepper,
      std::cref(vf), y, t_start, t_end, h0,
      observer);
  (void) steps; // TODO: do something with steps??

  State u = vf.decodeState(s.t()+dt, y);
  auto loads = vf.decodeLoads(s.t()+dt, y);
  return std::make_pair(u, loads);
}





/** Use an SDE integrator to evolve the State forward in time.
 */
std::pair<State, std::map<Transition*, double>>
    sde_integrate(double dt, const VectorField & vf, const State & s, double & h, Rng & rng) {
  // define the EM stepper
  EulerMaruyamaStepper stepper(rng); // rng is passed by ref
  // TODO: verify that s and vf are compatible
  RealVec y = vf.encode(s);

  double t_start = s.t();
  double t_end = s.t() + dt;
  double h0 = std::min(dt, h);
  size_t steps = odeint::integrate_const(stepper, std::cref(vf),
      y, t_start, t_end, h0);
  (void) steps; // TODO: do something with steps??

  State u = vf.decodeState(t_end, y);
  auto loads = vf.decodeLoads(t_end, y);
  return std::make_pair(u, loads);
}









Transition* choose_random_transition(const std::vector<Transition*> & transitions,
    const std::vector<double> & rates, double lambda, Rng & rng) {
  if ( lambda <= 0 ) {
    std::cerr << "# WARNING: lambda not positive: " << lambda << RIGHT_HERE << std::endl;
    return nullptr;
  }
  double u = 0.0;
  try {
    u = rng.uniform(0, lambda);
  } catch ( const std::exception & ex ) {
    std::cout << ex.what() << RIGHT_HERE << std::endl;
    throw ex;
  }
  auto transit = transitions.begin();
  for ( auto r : rates ) {
    u -= r;
    if ( u < 0 ) {
      break;
    } // else...
    ++transit;
  }
  if ( transit != transitions.end() ) {
    return *transit;
  } else {
    return nullptr; // signals failure
  }
}

Transition* choose_random_transition(const std::map<Transition*, double> & loads,
    double Lambda, Rng & rng) {
  if ( Lambda <= 0 ) {
    std::cerr << "# WARNING: lambda not positive: " << Lambda << RIGHT_HERE << std::endl;
    return nullptr;
  }
  double u = 0.0;
  try { // FIXME: Lambda can be nan!!
    u = rng.uniform(0, Lambda);
  } catch ( const std::exception & ex ) {
    std::cout << ex.what() << RIGHT_HERE << std::endl;
    throw ex;
  }
  //double u = rng.uniform(0, Lambda);
  for ( auto & load : loads ) {
    u -= load.second;
    if ( u < 0 ) {
      return load.first;
    }
  } // for load in loads
  return nullptr; // signals failure
}
