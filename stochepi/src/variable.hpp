#ifndef VARIABLE_HPP_
#define VARIABLE_HPP_

#include <iostream>

#include "aux.hpp"

/* The basic data type for the Hybrid model can be either
 * discrete or continuous
 */
class Variable : public Printable {
public:
  // default constructor: discrete 0
  Variable() : is_discrete(true), disc_value(0), cont_value(0.0) {}
  Variable(int x) : is_discrete(true), disc_value(x), cont_value(x) {}
  Variable(double x) : is_discrete(false), disc_value(x), cont_value(x) {}
  /// assinging with a double sets is_discrete to false
  const Variable & operator=(double x);
  /// assinging with an int sets is_discrete to true
  const Variable & operator=(int x);
  /// type conversion to double
  operator double() const;
  /// explicit value() method returns the xxxx_value cast to double
  double value() const;
  void make_discrete();
  void make_continuous();
  void print(std::ostream & os) const override;
  // data (TODO: make protected)
  bool is_discrete;
  int disc_value;
  double cont_value;
protected:
  /* empty */
};



#endif
