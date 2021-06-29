#ifndef DATA_HPP_
#define DATA_HPP_

/** @file data.hpp
 * @short Define timeseries data object and load data from files.
 * @todo refactor the load data functions
 */

#include <vector>
#include <iostream>
#include <string>

#include "aux.hpp"


class Observation : public Printable {
public:
  /** Constructor taking vectors for observations and censor codes.
   * The censor code vector is resized to the value vector,
   * and when it is shorter, UNCENSORED_CODE is added to the end
   */
  Observation(double t, std::vector<double> x, std::vector<int> c, std::string event);
  /** Constructor for when the observation is one dimensional
   */
  Observation(double t, double x, int c, std::string event) :
    t(t), x(std::vector<double>(1,x)), c(std::vector<int>(1,c)),
    event(event) { /* empty */ }
  double t; // time of observation
  std::vector<double> x; // the actual observation
  std::vector<int> c; // censor code
  std::string event;
  void print(std::ostream & os) const;
};


class Timeseries : public Printable {
public:
  std::vector<Observation> obs;
  void print(std::ostream & os) const;
  size_t size() const;
  std::string identifier;
};

/** remove data with observation time larger than tmax */
void pruneTimeseries(Timeseries & txs, double tmax);


Timeseries load_scalar_data(std::string filename);
Timeseries load_vector_data(std::string filename);

std::vector<Timeseries> load_scalar_panel_data(std::string filename);
std::vector<Timeseries> load_vector_panel_data(std::string filename);

std::pair<double, bool> get_event_time(const std::string & event, const Timeseries & txs);


#endif
