#include "parameter.hpp"

#include <algorithm> // for_each
#include <numeric> // accumulate
#include <sstream>

#include "macros.hpp"



Parameters::Parameters(int n) {
  paramvec.resize(n, 0.0);
  shape.resize(n,1);
}

Parameters::Parameters(const std::vector<size_t> & shape) : shape(shape) {
  size_t n = std::accumulate(shape.begin(), shape.end(), 0); // default accumulate is addition
  paramvec.resize(n, 0.0);
}



const ParPrior & Parameters::operator[](int i) const {
  // retrieve a parameter by it's index (can use enum)
  if ( i < 0 || i >= int(paramvec.size()) ) {
    std::stringstream message;
    message << "parameter index " << i << " out of range" << RIGHT_HERE;
    throw std::range_error(message.str());
  } // else ...
  return paramvec[i];
}

ParPrior & Parameters::operator[](int i) {
  // retrieve a parameter by it's index (can use enum)
  if ( i < 0 || i >= int(paramvec.size()) ) {
    std::stringstream message;
    message << "parameter index " << i << " out of range" << RIGHT_HERE;
    throw std::range_error(message.str());
  } // else ...
  return paramvec[i];
}

const ParPrior & Parameters::operator[](std::string name) const {
  for ( auto & param : paramvec ) {
    if ( param.getName() == name ) {
      return param;
    }
  } // else throw exception.
  std::stringstream message;
  message << "parameter with name " << name << " not found" << RIGHT_HERE;
  throw std::range_error(message.str());
}

ParPrior & Parameters::operator[](std::string name) {
  for ( auto & param : paramvec ) {
    if ( param.getName() == name ) {
      return param;
    }
  } // else throw exception.
  std::stringstream message;
  message << "parameter with name " << name << " not found" << RIGHT_HERE;
  throw std::range_error(message.str());
}





std::pair<size_t, size_t> Parameters::unflatten_index(size_t i) const {
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

size_t Parameters::flatten_index(size_t i, size_t j) const {
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


bool Parameters::isScalar(size_t i) const {
  if ( i >= shape.size() ) {
    throw std::range_error("index out of range" + RIGHT_HERE);
  } // else...
  return ( shape[i] == 1 );
}

const ParPrior & Parameters::operator()(size_t i, size_t j) const {
  return paramvec[flatten_index(i, j)];
}

ParPrior & Parameters::operator()(size_t i, size_t j) {
  return paramvec[flatten_index(i, j)];
}

// access to variables with single index can only be used for scalars
const ParPrior & Parameters::operator()(size_t i) const {
  if ( !isScalar(i) ) {
    throw std::logic_error("variable is not a scalar" + RIGHT_HERE);
  }
  return paramvec[flatten_index(i, 0)];
}

ParPrior & Parameters::operator()(size_t i) {
  if ( !isScalar(i) ) {
    throw std::logic_error("variable is not a scalar" + RIGHT_HERE);
  }
  return paramvec[flatten_index(i, 0)];
}



size_t Parameters::flat_size() const {
  return paramvec.size();
}

size_t Parameters::size() const {
  return shape.size();
}



void Parameters::mutate(Rng & rng, double rel_temp, double t) {
  // mutate standard parameters in the paramvec
  auto fun = [&](ParPrior & par) {
    if ( par.isRandomEffects() ) {
      // random walk with metropolis-hastings step: reject or accept
      double ll_old = par.loglike();
      par.mutate(rng, rel_temp, t);
      double u = 1 - rng.uniform(); // uniform(0, 1) includes 0, excludes 1
      if ( log(u) > par.loglike() - ll_old ) {
        par.reject();
      }
    } else {
      // simple random walk
      par.mutate(rng, rel_temp, t);
    }
  };
  std::for_each(paramvec.begin(), paramvec.end(), fun);
}

void Parameters::select(int r) {
  auto fun = [r](ParPrior & par) {
    if ( par.size() > 1 ) par.select(r);
  };
  std::for_each(paramvec.begin(), paramvec.end(), fun);
}

void Parameters::deselect() {
  auto fun = [](ParPrior & par){par.deselect();};
  std::for_each(paramvec.begin(), paramvec.end(), fun);
}

void Parameters::lockAllBut(int r) {
  // locks all parameter elements except the r-th if at least one element is NOT locked
  auto fun = [r](ParPrior & par) {
    if ( par.size() > 1 && !par.isLocked() ) par.lockAllBut(r);
  };
  std::for_each(paramvec.begin(), paramvec.end(), fun);
}

void Parameters::removeSingleLocks() {
  /* undoes the effect of lockAllBut. That is:
   * removes all locks for random effects parameters
   * if at least one element is NOT locked
   */
  auto fun = [](ParPrior & par){par.removeSingleLocks();};
  std::for_each(paramvec.begin(), paramvec.end(), fun);
}


double Parameters::loglike() const { // prior likelihood
  double ll = 0.0;
  for ( auto & par : paramvec ) {
    ll += par.loglike();
  }
  return ll;
}


void Parameters::print(std::ostream & os) const {
  // use lambda to avoid repeating code...
  auto fun = [&](const ParPrior & par) {
    if ( !par.isLocked() ) {
      os << par << std::endl;
      if ( par.isRandomEffects() ) {
        os << *(par.getLoc()) << std::endl;
        os << *(par.getScale()) << std::endl;
      }
    }
  };

  os << "<parameters >" << std::endl;
  std::for_each(paramvec.begin(), paramvec.end(), fun);
  os << "</parameters>";
}
