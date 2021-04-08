#include "model.hpp"

#include <sstream>



Model::Model(const State & s0, const std::vector<Transition*> transitions,
    const Parameters & par, bool sde, const std::string & name) :
    s0(s0), par(par), sde(sde), name(name) {
  replace_transitions(transitions);
  // defaults for functions...
  likfun = [](const Observation & obs, const State & state,
      const Parameters & par) -> std::pair<double, bool> {
    return std::make_pair(0.0, true);
  };
  initfun = [s0](const Parameters & par, Rng & rng) -> State {
    return s0;
  };
  coftestfun = [](const State & state) -> bool {
    return false;
  };
  coflikfun = [](const Observation & obs) -> std::pair<double, bool> {
    return std::make_pair(0.0, true);
  };
  obsfun = [](const State & state, const Parameters & par) -> State {
    State obs;
    obs.t() = state.t();
    return obs; // no observables by default: empty state object
  };
  modfun = [](State & state, const Observation & obs, const Parameters & par,
      Rng & rng) {
    // be default: don't modify the state...
  };
}

Model::Model(const Model & model) : s0(model.s0), par(model.par),
  likfun(model.likfun), initfun(model.initfun), coftestfun(model.coftestfun),
  coflikfun(model.coflikfun), obsfun(model.obsfun), modfun(model.modfun),
  sde(model.sde) {
  replace_transitions(model.transitions);
}

Model & Model::operator=(const Model & model) {
  if ( this != &model ) {
    s0 = model.s0;
    par = model.par;
    likfun = model.likfun;
    initfun = model.initfun;
    coftestfun = model.coftestfun;
    coflikfun = model.coflikfun;
    obsfun = model.obsfun;
    modfun = model.modfun;
    sde = model.sde;
    name = model.name;
    replace_transitions(model.transitions);
  }
  return *this;
}


/** Make a representation of the model in the dot language.
 * @todo: make more informative representations with variables as nodes
 */
std::string Model::draw() const {
  std::stringstream ss_graph;
  // build the label for the center node
  std::stringstream ss_state;
  std::string sep = "";
  for ( size_t i = 0; i < s0.size(); ++i ) {
    ss_state << sep << s0.getName(i);
    sep = ", ";
  }

  ss_graph << "digraph model {" << std::endl;
  ss_graph << "  s0 [label=\"" << ss_state.str() << "\"]" << std::endl;

  for ( size_t i = 0; i < transitions.size(); ++i ) {
    State sprime = transitions[i]->apply(s0);
    std::stringstream ss_state;
    std::string sep = "";
    for ( size_t j = 0; j < s0.size(); ++j ) {
      int incr = sprime(j, 0).disc_value - s0(j, 0).disc_value;
      // FIXME: the 0 is a quick fix... properly handle vector-valued variables
      std::string addl;
      if ( incr > 0 ) {
        addl = "+" + std::to_string(incr);
      } else if ( incr < 0 ) {
        addl = "-" + std::to_string(-incr);
      } // else zero: no addl
      ss_state << sep << s0.getName(j) << addl;
      sep = ", ";
    }
    ss_graph << "  s" << i+1 << " [label=\"" << ss_state.str() << "\"]" << std::endl;
  }


  for ( size_t i = 0; i < transitions.size(); ++i ) {
    ss_graph << "  s0 -> s" << i+1 << " [label=\"" << transitions[i]->getName() << "\"]" << std::endl;
  }

  ss_graph << "}";
  return ss_graph.str();
}


std::string Model::getName() const {
  return name;
}

const std::vector<Transition*> & Model::getTransitions() const {
  return transitions;
}
