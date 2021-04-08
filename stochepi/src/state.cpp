#include "state.hpp"

#include <algorithm>
#include <numeric>
#include <sstream>

#include "macros.hpp"


State::State(size_t n) : _t(0.0) {
  vars.resize(n, 0);
  shape.resize(n, 1);
  names.resize(n); // TODO: good default names?
}

State::State(const std::vector<std::string> & names) : _t(0.0), names(names) {
  vars.resize(names.size(), 0);
  shape.resize(names.size(), 1);
}

State::State(const std::vector<std::string> & names,
    const std::vector<size_t> & shape) : _t(0.0), shape(shape), names(names) {
  // some checks
  if ( shape.size() != names.size() ) {
    throw std::invalid_argument("names and shape must have the same size" + RIGHT_HERE);
  } // else...
  size_t n = std::accumulate(shape.begin(), shape.end(), 0); // default accumulate is addition
  vars.resize(n, 0);
}


State::State(const std::vector<std::string> & names, size_t shape) : _t(0.0),
    names(names) {
  this->shape.assign(names.size(), shape);
  size_t n = shape * names.size();
  // default accumulate is addition
  vars.resize(n, 0);
}



State State::zero_like() const { // returns a zero with the same shape (and names)
  return State(names, shape);
}

std::vector<double> State::vec() const {
  std::vector<double> y(vars.size());
  for ( size_t i = 0; i < vars.size(); ++i ) {
    y[i] = vars[i].value();
  }
  return y;
}



const double & State::t() const {
  return _t;
}

double & State::t() {
  return _t;
}



size_t State::size() const {
  return shape.size();
}

size_t State::size(size_t i) const { // returns shape[i]
  if ( i >= shape.size() ) {
    throw std::range_error("index out of range" + RIGHT_HERE);
  } // else...
  return shape[i];
}

size_t State::flat_size() const {
  return vars.size();
}

bool State::empty() const {
  return shape.empty();
}

const std::string & State::getName(size_t i) const {
  if ( i >= names.size() ) {
    throw std::range_error("index out of range" + RIGHT_HERE);
  } // else...
  return names[i];
}



const Variable & State::operator[](size_t i) const {
  if ( i >= vars.size() ) {
    throw std::range_error("index out of range" + RIGHT_HERE);
  } // else...
  return vars[i];
}

Variable & State::operator[](size_t i) {
  if ( i >= vars.size() ) {
    throw std::range_error("index out of range" + RIGHT_HERE);
  } // else...
  return vars[i];
}



bool State::isNonNegative() const {
  // test if discrete variables are non-negative
  for ( auto & x : vars ) {
    if ( x.is_discrete && x.disc_value < 0 ) {
      return false;
    }
  } // loop over variables
  return true;
}



bool State::compareShapes(const State & s) const {
  if ( flat_size() != s.flat_size() || size() != s.size() ) {
    return false;
  }
  for ( size_t i = 0; i < shape.size(); ++i ) {
    if ( shape[i] != s.shape[i] ) return false;
  }
  return true;
}


bool State::compareSignatures(const State & s) const {
  /* returns true if all the variables of this and s have the same type
   * and false otherwise. The flat_size()s must agree as well.
   */
  if ( vars.size() != s.vars.size() ) {
    return false;
  } // else...
  for ( size_t i = 0; i < vars.size(); ++i ) {
    if ( vars[i].is_discrete != s.vars[i].is_discrete ) {
      return false;
    }
  } // for loop over vars and s.vars
  // if this point is reached, all variables have the same type and shapes agree
  return true;
}


void State::print(std::ostream & os) const {
  size_t k = 0; // index used for vars (flat index)
  os << "<state t='" << _t << "' >" << std::endl;
  for ( size_t i = 0; i < shape.size(); ++i ) {
    os << "<var_vec name='" << names[i] << "' "
       << "dim='" << shape[i] << "' "
       << ">" << std::endl;
    for ( size_t j = 0; j < shape[i]; ++j ) {
      os << vars[k] << std::endl;
      ++k; // inner loop: increment flat index k each time
    }
    os << "</var_vec>" << std::endl;
  }
  os << "</state>";
}


std::string State::rdump() const { // represent the state in the R dump format
  /* the type of the variables (integer, real,  ...)
   * is coded as follows:
   * int: 0, real: 1
   */
  std::stringstream ss;
  size_t k = 0; // index used for vars (flat index)
  ss << "t <- " << _t << std::endl;
  for ( size_t i = 0; i < shape.size(); ++i ) {
    if ( shape[i] == 1 ) {
      ss << names[i] << " <- " << vars[k].value() << std::endl;
      ss << names[i] << ".type <- " << (vars[k].is_discrete ? 0 : 1 ) << std::endl;
      ++k; // inner loop: increment flat index k each time
    } else {
      // TODO: give types of vector-valued variables
      ss << names[i] << " <- c(";
      std::string sep = "";
      for ( size_t j = 0; j < shape[i]; ++j ) {
        ss << sep << vars[k].value() << std::endl;
        ++k; // inner loop: increment flat index k each time
        if ( sep.empty() ) sep = ", ";
      }
      ss << ")" << std::endl;
    }
  }
  return ss.str();
}



std::pair<size_t, size_t> State::unflatten_index(size_t i) const {
  for ( size_t j = 0; j < shape.size(); ++j ) {
    if ( i < shape[j] ) {
      return std::make_pair(j, i);
    } else {
      i -= shape[j];
    }
  }
  // if this point can be reached, the index is out of range
  throw std::range_error("given index too large" + RIGHT_HERE);
}

size_t State::flatten_index(size_t i, size_t j) const {
  if ( i >= shape.size() ) {
    throw std::range_error("first index out of range" + RIGHT_HERE);
  } // else...
  if ( j >= shape[i] ) {
    throw std::range_error("second index out of range" + RIGHT_HERE);
  } // else...
  int k = 0;
  for ( size_t ii = 0; ii < i; ++ii ) {
    k += shape[ii];
  }
  return k + j;
}


bool State::isScalar(size_t i) const {
  if ( i >= shape.size() ) {
    throw std::range_error("index out of range" + RIGHT_HERE);
  } // else...
  return ( shape[i] == 1 );
}

const Variable & State::operator()(size_t i, size_t j) const {
  return vars[flatten_index(i, j)];
}

Variable & State::operator()(size_t i, size_t j) {
  return vars[flatten_index(i, j)];
}

// access to variables with single index can only be used for scalars
const Variable & State::operator()(size_t i) const {
  if ( !isScalar(i) ) {
    throw std::logic_error("variable is not a scalar" + RIGHT_HERE);
  }
  return vars[flatten_index(i, 0)];
}

Variable & State::operator()(size_t i) {
  if ( !isScalar(i) ) {
    throw std::logic_error("variable is not a scalar" + RIGHT_HERE);
  }
  return vars[flatten_index(i, 0)];
}



State::iterator State::begin() {
  return vars.begin();
}

State::iterator State::end() {
  return vars.end();
}

State::const_iterator State::begin() const {
  return vars.begin();
}

State::const_iterator State::end() const {
  return vars.end();
}



State State::concat(const State & s) const {
  State con(*this); // copy of this
  // add data from s
  if ( !s.empty() ) {
    con.vars.insert(con.vars.end(), s.vars.begin(), s.vars.end());
    con.shape.insert(con.shape.end(), s.shape.begin(), s.shape.end());
    con.names.insert(con.names.end(), s.names.begin(), s.names.end());
  }
  return con;
}

State concat(const State & s1, const State & s2) {
  return s1.concat(s2);
}




bool compareSignatures(const State & s1, const State & s2) {
  return s1.compareSignatures(s2);
}

bool compareShapes(const State & s1, const State & s2) {
  return s1.compareShapes(s2);
}
