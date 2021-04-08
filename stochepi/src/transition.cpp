#include "transition.hpp"

#include <algorithm> // for_each


Transition* Transition::dup() const {
  return new Transition(*this); // SHOULD THIS BE REDEFINED FOR CHILDREN
}
State Transition::apply(const State & s) const {
  return s; // do nothing
}
double Transition::rate(const State & , const Parameters & ) const {
  return 0.0; // nothing never happens
}
std::string Transition::getName() const {
  return name;
}
void Transition::print(std::ostream & os) const {
  os << getName();
}
double & Transition::load() {
  return _load;
}
double Transition::load() const {
  return _load;
}


// methods for TransitionList

TransitionList::~TransitionList() {
  clear_transitions();
}

void TransitionList::replace_transitions(const std::vector<Transition*> & new_transitions) {
  // delete old Transitions
  clear_transitions();
  // resize old vector
  transitions.reserve(new_transitions.size());
  // duplicate the transitions and put them in the vector
  for ( auto trans : new_transitions ) {
    transitions.push_back(trans->dup());
  }
}

void TransitionList::clear_transitions() {
  std::for_each(transitions.begin(), transitions.end(),
      [](Transition* trans){delete trans;});
  transitions.clear();
}
