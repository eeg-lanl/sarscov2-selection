#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <string>
#include <functional>

#include "state.hpp"
#include "transition.hpp"
#include "data.hpp"


/** Function that returns the log-likelihood of an observation,
 * given a state and  parameters.
 */
typedef std::function<std::pair<double, bool>(const Observation&,
    const State&, const Parameters&)> LikelihoodFunction;


/** Generate an initial state using an RNG and parameters
 */
typedef std::function<State(const Parameters&, Rng&)>
    InitialStateFunction;

/** Generate observables from the State and Parameters
 */
typedef std::function<State(const State&, const Parameters&)>
    ObservableFunction;


/** A model can define a "coffin" state (e.g. extinct virus population)
 * That we have to keep track of, but does not evolve any further.
 */
typedef std::function<bool(const State&)> CoffinTestFunction;


/** for the coffin state, we have to compute a likelihood
 * of an observation */
typedef std::function<std::pair<double, bool>(const Observation&)>
    CoffinLikFunction;


/** Allow for events to change the state of the model at observation times.
 * e.g. reset a cumulative value used for the next observation.
 */
typedef std::function<void(State&, const Observation&, const Parameters&, Rng&)>
    StateModifier;


/** A Model contains:
 * 1) a list of Transitions (inherited from TransitionList)
 * 2) an initial State,
 * 3) default Parameters
 * 3) a likelihood function for observations
*/
class Model : public TransitionList {
public:
  Model() { /* empty */ }
  Model(const State & s0, const std::vector<Transition*> transitions,
    const Parameters & par, bool sde, const std::string & name);
  Model(const Model & model);
  Model & operator=(const Model & model);
  virtual ~Model() { /* empty */ };
  State s0;
  Parameters par;
  LikelihoodFunction likfun;
  InitialStateFunction initfun;
  CoffinTestFunction coftestfun;
  CoffinLikFunction coflikfun;
  ObservableFunction obsfun;
  StateModifier modfun;
  bool sde; // determines kind of model (TODO: expand options)
  // make a graphical representation of the model
  std::string draw() const;
  std::string getName() const;
  const std::vector<Transition*> & getTransitions() const;
protected:
  std::string name;
  // addl data (or make data protected)
};


#endif
