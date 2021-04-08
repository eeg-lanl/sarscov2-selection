#ifndef PARSE_PARAMS_HPP
#define PARSE_PARAMS_HPP

/** @file parse_params.hpp
 * @short declare functions to parse parameter files
 */

#include <list>
#include <string>
#include <iostream>


struct ParSpecs {
  std::string name;
  double value;
  bool is_locked;
  double sigma_rw;
  bool is_lbound;
  double lbound;
  bool is_ubound;
  double ubound;
};

std::ostream & operator<<(std::ostream & os, const ParSpecs & parSpecs);

/** Parse the content of a parameter file
 *
 * The columns of the parameter file should be as the following examples:
 *
 * # Parameter between 0 and 2.5 with random walk sd 0.1:
 * beta 0.1 0.0 2.5
 * # Unbounded parameter with random walk sd 1.0
 * t0 1.0 -inf inf
 * # Positive parameter that is locked
 * gamma 0 0.0 inf
 *
 */
std::list<ParSpecs> parse_param_file(std::string filename);


#endif
