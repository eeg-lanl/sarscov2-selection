#include "parse_params.hpp"

#include <fstream>
#include <sstream>
#include <exception>
#include <algorithm>

#include "macros.hpp"


std::ostream & operator<<(std::ostream & os, const ParSpecs & psc) {
  os << "<par_spec "
     << "name='" << psc.name << "' "
     << "value='" << psc.value << "' "
     << "sigma_rw='" << psc.sigma_rw << "' "
     << "is_locked='" << std::boolalpha << psc.is_locked << std::noboolalpha << "' "
     << "lbound='" << psc.lbound << "' "
     << "is_lbound='" << std::boolalpha << psc.is_lbound << std::noboolalpha << "' "
     << "ubound='" << psc.ubound << "' "
     << "is_ubound='" << std::boolalpha << psc.is_ubound << std::noboolalpha << "' "
     << "/>";
  return os;
}


std::list<ParSpecs> parse_param_file(std::string filename) {
  std::ifstream file(filename);
  if ( !file.good() ) {
    throw std::runtime_error("unable to open file " + filename + RIGHT_HERE);
  }
  std::list<ParSpecs> parSpecsList;
  std::string line;
  while ( std::getline(file, line) ) {
    if ( line.empty() || line[0] == '#' ) continue;
    // else: line is not empty or a comment
    // get parameter name
    auto start = find_if_not(line.begin(), line.end(), isspace);
    auto stop = find_if(start, line.end(), isspace);
    std::string parname(start, stop);
    // get value
    start = find_if_not(stop, line.end(), isspace);
    stop = find_if(start, line.end(), isspace);
    std::string value(start, stop);
    // get sigma_rw
    start = find_if_not(stop, line.end(), isspace);
    stop = find_if(start, line.end(), isspace);
    std::string sigma_rw(start, stop);
    // get lower bound
    start = find_if_not(stop, line.end(), isspace);
    stop = find_if(start, line.end(), isspace);
    std::string lbound(start, stop);
    // get upper bound
    start = find_if_not(stop, line.end(), isspace);
    stop = find_if(start, line.end(), isspace);
    std::string ubound(start, stop);
    // create a ParSpecs object
    ParSpecs parSpecs;
    parSpecs.name = parname;
    // use stringstreams to convert strings to values
    std::stringstream(value) >> parSpecs.value;
    std::stringstream(sigma_rw) >> parSpecs.sigma_rw;
    // sigma_rw <= 0 means that the parameter is fixed
    parSpecs.is_locked = (parSpecs.sigma_rw <= 0);
    // lower bound
    if ( lbound == "-inf" ) {
      parSpecs.is_lbound = false;
      parSpecs.lbound = 0.0; // dummy value
    } else {
      parSpecs.is_lbound = true;
      std::stringstream(lbound) >> parSpecs.lbound;
    }
    // upper bound
    if ( ubound == "inf" ) {
      parSpecs.is_ubound = false;
      parSpecs.ubound = 0.0; // dummy value
    } else {
      parSpecs.is_ubound = true;
      std::stringstream(ubound) >> parSpecs.ubound;
    }
    // add specs to the map
    parSpecsList.push_back(parSpecs);
  }
  return parSpecsList;
}
