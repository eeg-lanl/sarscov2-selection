#ifndef STATE_HPP_
#define STATE_HPP_

#include <vector>
#include <string>

#include "aux.hpp"
#include "variable.hpp"

/** @todo named variables (pass a list of strings at construction)
 * overload operator[] to get access to variables using their names
 */
class State : public Printable {
public:
  State() : _t(0.0) { /* 0-dimensional State */ }
  State(size_t ); // n-dimensional state
  template<class Numeric>
  State(double t, const std::vector<Numeric> & y) : _t(t) {
    shape.resize(y.size(), 1);
    vars.resize(y.size());
    for ( size_t i = 0; i < vars.size(); ++i ) {
      vars[i] = y[i];
    }
  }
  State(const std::vector<std::string> & names); // TODO variadic constructor??
  State(const std::vector<std::string> & names, const std::vector<size_t> & shape);
  /// Repeat shape for all vars
  State(const std::vector<std::string> & names, size_t shape);
  State zero_like() const; // returns a zero with the same shape (and names)
  std::vector<double> vec() const;
  size_t size() const; // returns shape.size()
  size_t size(size_t i) const; // returns shape[i]
  size_t flat_size() const; // returns vars.size() = sum_i shape[i]
  bool empty() const; // size = 0
  const std::string & getName(size_t i) const;
  void print(std::ostream & os) const override;
  std::string rdump() const; // represent the state in the R dump format
  bool isNonNegative() const;
  // check that the States have the same shape
  bool compareShapes(const State & ) const;
  // check that the variables of this and an other state have identical type
  bool compareSignatures(const State & ) const;
  bool isScalar(size_t i) const; // check that the i-th block is a scalar
  const Variable & operator[](size_t i) const; // access variable directly
  Variable & operator[](size_t i);
  const Variable & operator()(size_t i, size_t j) const;
  Variable & operator()(size_t i, size_t j); // access elements according to shape
  // access to variables with single index can only be used for scalars
  const Variable & operator()(size_t i) const;
  Variable & operator()(size_t i); // access elements according to shape
  typedef std::vector<Variable>::iterator iterator;
  typedef std::vector<Variable>::const_iterator const_iterator;
  iterator begin();
  const_iterator begin() const;
  iterator end();
  const_iterator end() const;
  const double & t() const;
  double & t();
  // concatentate states
  State concat(const State & s) const;
protected:
  std::pair<size_t, size_t> unflatten_index(size_t i) const;
  size_t flatten_index(size_t i, size_t  j) const;
  // data
  double _t; // time: TODO: consider removing time from State
  std::vector<Variable> vars; // e.g. A, I, V
  std::vector<size_t> shape;
  std::vector<std::string> names;
};

// check that the States have the same shape
bool compareShapes(const State & s1, const State & s2);

// compareSIgnatures calls the State method with the same name
bool compareSignatures(const State & s1, const State & s2);

// external concatenation function
State concat(const State & s1, const State & s2);



#endif
