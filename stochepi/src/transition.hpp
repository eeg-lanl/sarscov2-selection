#ifndef TRANSITION_HPP_
#define TRANSITION_HPP_

#include <string>
#include <vector>
#include <map> // for parameter name maps

#include "aux.hpp"
#include "state.hpp"
#include "parameter.hpp"


class Transition : public Indexable, public Printable {
public:
  Transition() : name("Identity"), _load(0.0) { /* empty */ }
  Transition(std::string name) : name(name), _load(0.0) { /* empty */ }
  virtual ~Transition() { /* empty */ }
  virtual Transition* dup() const;
  virtual State apply(const State & ) const; // TODO: in-place mutation of state
  virtual double rate(const State & , const Parameters & ) const;
  std::string getName() const;
  void print(std::ostream & ) const;
  double & load();
  double load() const;
protected:
  std::string name;
private:
  double _load;
};

/** Auxiliary class to give classes derived from Transition the correct dup method
 * FIXME: mimick covariant return types of the dub method with some auxiliary
 * function
 * USAGE:
 * class DerivedTransition : public CloneableTransition<DerivedTransition> {
 *   DerivedTransition() : CloneableTransition<DerivedTransition>("NameOfDerivedTransition") { ... }
 *   // implement apply and rate
 * };
 */
template<class Derived>
class CloneableTransition : public Transition {
public:
  CloneableTransition(std::string name) : Transition(name) { /* empty */ }
  virtual ~CloneableTransition() { /* empty */ };
  virtual Transition* dup() const override {
    return new Derived(*static_cast<const Derived*>(this));
  }
};


/* a class wrapper for std::vector<Transition*> that takes care of
 * correctly making deep copies of the Transitions.
 * TODO: make a constructor with an initializer list.
 * Is this safe? e.g.
 */
class TransitionList {
public:
  virtual ~TransitionList();
protected:
  // aux function for constructor and copy constructor
  void replace_transitions(const std::vector<Transition*> & );
  void clear_transitions();
  // member
  std::vector<Transition*> transitions;
};



#endif
