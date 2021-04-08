#include "data.hpp"

#include <fstream>
#include <sstream>
#include <exception>
#include <algorithm>
#include <list>

#include "macros.hpp"

Observation::Observation(double t, std::vector<double> x, std::vector<int> c,
    std::string event) :  t(t), x(x), c(c), event(event) {
  // make sure that x and c have the same size, adding UNCENSORED_CODE where needed
  c.resize(x.size(), UNCENSORED_CODE);
}


void Observation::print(std::ostream & os) const {
  if ( x.size() != c.size() ) {
    throw std::logic_error("invalid observation" + RIGHT_HERE);
  }
  os << t;
  for ( size_t i = 0; i < x.size(); ++i ) {
    os << "\t" << x[i] << "\t" << c[i];
  }
  os << "\t" << event;
}

void Timeseries::print(std::ostream & os) const {
  for ( auto & ob : obs ) {
    os << identifier << "\t" << ob << std::endl;
  }
}

size_t Timeseries::size() const {
  return obs.size();
}

Timeseries load_scaler_data(std::string filename) {
  Timeseries txs;
  std::ifstream file(filename.c_str());
  if ( !file.is_open() ) {
    throw std::runtime_error("cannot open file '" + filename + "'" + RIGHT_HERE);
  } // else...
  std::string line;
  while ( std::getline(file, line) ) {
    if ( line.empty() || line[0] == '#' ) continue;
    // else: line is not empty or a comment
    // get ID
    auto start = find_if_not(line.begin(), line.end(), isspace);
    auto stop = find_if(start, line.end(), isspace);
    std::string identifier(start, stop);
    // get time point
    start = find_if_not(stop, line.end(), isspace);
    stop = find_if(start, line.end(), isspace);
    std::string time(start, stop);
    // get value
    start = find_if_not(stop, line.end(), isspace);
    stop = find_if(start, line.end(), isspace);
    std::string val(start, stop);
    // get censoring type
    start = find_if_not(stop, line.end(), isspace);
    stop = find_if(start, line.end(), isspace);
    std::string cens(start, stop);
    // get event
    start = find_if_not(stop, line.end(), isspace);
    stop = find_if(start, line.end(), isspace);
    std::string event(start, stop);
    // use stringstreams to convert strings to values
    double t;
    std::stringstream(time) >> t;
    double x;
    std::stringstream(val) >> x;
    int c;
    std::stringstream(cens) >> c;

    Observation ob(t, x, c, event);
    txs.obs.push_back(ob);
    // add identifier to the timeseries
    txs.identifier = identifier;
  }
  return txs;
}

Timeseries load_vector_data(std::string filename) {
  Timeseries txs;
  std::ifstream file(filename.c_str());
  if ( !file.is_open() ) {
    throw std::runtime_error("cannot open file '" + filename + "'" + RIGHT_HERE);
  } // else...
  std::string line;
  while ( std::getline(file, line) ) {
    if ( line.empty() || line[0] == '#' ) continue;
    // else: line is not empty or a comment
    // get ID
    auto start = find_if_not(line.begin(), line.end(), isspace);
    auto stop = find_if(start, line.end(), isspace);
    std::string identifier(start, stop);
    // get time point
    start = find_if_not(stop, line.end(), isspace);
    stop = find_if(start, line.end(), isspace);
    std::string time(start, stop);
    // get event
    start = find_if_not(stop, line.end(), isspace);
    stop = find_if(start, line.end(), isspace);
    std::string event(start, stop);
    // use stringstreams to convert strings to values
    double t;
    std::stringstream(time) >> t;
    // get values and censor code
    std::vector<double> xs;
    std::vector<int> cs;
    while ( stop != line.end() ) {
      start = find_if_not(stop, line.end(), isspace);
      stop = find_if(start, line.end(), isspace);
      std::string val(start, stop);
      double x;
      std::stringstream(val) >> x;
      xs.push_back(x);
      // get censoring type
      start = find_if_not(stop, line.end(), isspace);
      stop = find_if(start, line.end(), isspace);
      std::string cens(start, stop);
      int c;
      std::stringstream(cens) >> c;
      cs.push_back(c);
    }
    // create a new observation
    Observation ob(t, xs, cs, event);
    txs.obs.push_back(ob);
    // add identifier to the timeseries
    txs.identifier = identifier;
  }
  return txs;
}




std::vector<Timeseries> load_scalar_panel_data(std::string filename) {
  std::list<Timeseries> txss;
  Timeseries txs; // the first entry
  std::ifstream file(filename.c_str());
  if ( !file.is_open() ) {
    throw std::runtime_error("cannot open file '" + filename + "'" + RIGHT_HERE);
  } // else...
  std::string line;
  while ( std::getline(file, line) ) {
    if ( line.empty() || line[0] == '#' ) continue;
    // else: line is not empty or a comment
    // get ID
    auto start = find_if_not(line.begin(), line.end(), isspace);
    auto stop = find_if(start, line.end(), isspace);
    std::string identifier(start, stop);
    // get time point
    start = find_if_not(stop, line.end(), isspace);
    stop = find_if(start, line.end(), isspace);
    std::string time(start, stop);
    // get value
    start = find_if_not(stop, line.end(), isspace);
    stop = find_if(start, line.end(), isspace);
    std::string val(start, stop);
    // get censoring type
    start = find_if_not(stop, line.end(), isspace);
    stop = find_if(start, line.end(), isspace);
    std::string cens(start, stop);
    // get event
    start = find_if_not(stop, line.end(), isspace);
    stop = find_if(start, line.end(), isspace);
    std::string event(start, stop);
    // use stringstreams to convert strings to values
    double t;
    std::stringstream(time) >> t;
    double x;
    std::stringstream(val) >> x;
    int c;
    std::stringstream(cens) >> c;
    // make new observation
    Observation ob(t, x, c, event);
    // add the observation to current or next timeseries
    if ( txs.identifier == identifier ) {
      // add obs to current timeseries
      txs.obs.push_back(ob);
    } else {
      if ( !txs.obs.empty() ) { // don't add the first empty timeseries
        txss.push_back(std::move(txs));
      }
      txs = Timeseries(); // start new empty timeseries
      txs.obs.push_back(std::move(ob));
      txs.identifier = identifier;
    }
  }
  // add the final timeseries to txss
  txss.push_back(std::move(txs));
  return std::vector<Timeseries>(txss.begin(), txss.end());
}



std::vector<Timeseries> load_vector_panel_data(std::string filename) {
  std::list<Timeseries> txss;
  Timeseries txs; // the first entry
  std::ifstream file(filename.c_str());
  if ( !file.is_open() ) {
    throw std::runtime_error("cannot open file '" + filename + "'" + RIGHT_HERE);
  } // else...
  std::string line;
  while ( std::getline(file, line) ) {
    if ( line.empty() || line[0] == '#' ) continue;
    // else: line is not empty or a comment
    // get ID
    auto start = find_if_not(line.begin(), line.end(), isspace);
    auto stop = find_if(start, line.end(), isspace);
    std::string identifier(start, stop);
    // get time point
    start = find_if_not(stop, line.end(), isspace);
    stop = find_if(start, line.end(), isspace);
    std::string time(start, stop);
    // get event
    start = find_if_not(stop, line.end(), isspace);
    stop = find_if(start, line.end(), isspace);
    std::string event(start, stop);
    // use stringstreams to convert strings to values
    double t;
    std::stringstream(time) >> t;
    // get values and censor code
    std::vector<double> xs;
    std::vector<int> cs;
    while ( stop != line.end() ) {
      start = find_if_not(stop, line.end(), isspace);
      stop = find_if(start, line.end(), isspace);
      std::string val(start, stop);
      double x;
      std::stringstream(val) >> x;
      xs.push_back(x);
      // get censoring type
      start = find_if_not(stop, line.end(), isspace);
      stop = find_if(start, line.end(), isspace);
      std::string cens(start, stop);
      int c;
      std::stringstream(cens) >> c;
      cs.push_back(c);
    }
    // create a new observation
    Observation ob(t, xs, cs, event);
    // add the observation to current or next timeseries
    if ( txs.identifier == identifier ) {
      // add obs to current timeseries
      txs.obs.push_back(ob);
    } else {
      if ( !txs.obs.empty() ) { // don't add the first empty timeseries
        txss.push_back(std::move(txs));
      }
      txs = Timeseries(); // start new empty timeseries
      txs.obs.push_back(std::move(ob));
      txs.identifier = identifier;
    }
  }
  // add the final timeseries to txss
  txss.push_back(std::move(txs));
  return std::vector<Timeseries>(txss.begin(), txss.end());
}

/* auxiliary function to get event time from timeseries.
 * This can fail whenever the event is not in the Timeseries.
 */
std::pair<double, bool> get_event_time(const std::string & event, const Timeseries & txs) {
  for ( auto obs : txs.obs ) {
    if ( obs.event == event ) {
      return std::make_pair(obs.t, true);
    }
  }
  return std::make_pair(0.0, false); // failure
}
