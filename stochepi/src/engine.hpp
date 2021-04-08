#ifndef ENGINE_HPP_
#define ENGINE_HPP_

#include <map>
#include <vector>
#include <functional>

#include "rng.hpp"
#include "aux.hpp"
#include "state.hpp"
#include "transition.hpp"
#include "vectorfield.hpp"
#include "macros.hpp"


constexpr double THRESHOLD = 50; // TODO: parameter
constexpr double WINDOW = 5; // TODO: parameter

/** Function object that determines the type of the state variables:
 * Either discrete or continuous
 */
typedef std::function<std::vector<bool>(const State&, const Parameters&)> VarTypeSwitcher;


class Engine : public TransitionList {
public:
  Engine(const State & , const std::vector<Transition*> , const Parameters & );
  Engine(const Engine & ); // copy constructor: duplicate transitions
  Engine & operator=(const Engine & ); // copy assignment constructor
  void reset(const State & s0, const Parameters & par);
  // evolve functions are not const as they modify h
  State evolve_gillespie(State s, double tmax, const Parameters & p, Rng & rng);
  State evolve_tauleap(State s, double tmax, const Parameters & p, Rng & rng);
  State evolve_ode(State s, double tmax, const Parameters & p);
  State evolve_sde(State s, double tmax, const Parameters & p, Rng & rng);
  /** The implementation of the "hybrid" model.
   *
   * Use an adaptive tau-leaping algorithm to evolve small populations,
   * and switch to ODEs (SDEs) for larger populations.
   */
  State evolve_hybrid(State s, double tmax, const Parameters & p,
      Rng & rng, bool sde=false);
  State evolve_hybrid(State s, double tmax, const Parameters & p,
      const VarTypeSwitcher & switcher, Rng & rng, bool sde=false);
  // static members
  const static double dt_max, dt_init, h_init;
protected:
  State s0;
  double h; // integration step
  double dt; // time step for tau-leaping
  VectorField vf;
};

std::pair<State, std::map<Transition*, double>>
    integrate(double , const VectorField & , const State & , double & );


std::pair<State, std::map<Transition*, double>>
    sde_integrate(double , const VectorField & , const State & , double & , Rng & );


// used in evolve_hybrid
Transition* choose_random_transition(const std::map<Transition*, double> & loads,
    double Lambda, Rng & rng);

// used in evolve_gillespie and evolve_tauleap
Transition* choose_random_transition(const std::vector<Transition*> & ,
    const std::vector<double> & , double , Rng & );


#endif
