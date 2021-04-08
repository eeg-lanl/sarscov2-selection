#include "splines.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

#include "macros.hpp"

using namespace splines;

size_t splines::num_knots(size_t order, size_t num_break_points) {
  return num_break_points + 2*(order-1);
}

size_t splines::num_basis(size_t order, size_t num_break_points) {
  return num_break_points + order - 2;
}

void splines::gen_knots_ref(size_t order, const std::vector<double> & break_points,
      std::vector<double> & knots) {
  size_t deg = order-1;
  size_t K = num_knots(order, break_points.size());
  knots.resize(K);
  // repeat end points
  for ( size_t j = 0; j < order-1; ++j ) {
    knots[j] = break_points.front();
    knots[K-1-j] = break_points.back();
  }
  // copy middle
  std::copy(break_points.begin(), break_points.end(), knots.begin()+deg);
}

void splines::gen_unif_knots_ref(size_t order, size_t num_break_points, double a, double b,
      std::vector<double> & knots) {
  size_t deg = order-1;
  size_t K = num_break_points + 2*deg;
  knots.resize(K);
  // repeat end points
  for ( size_t j = 0; j < order-1; ++j ) {
    knots[j] = a;
    knots[K-1-j] = b;
  }
  double h = (b-a) / (num_break_points - 1);
  for ( size_t j = 0; j < num_break_points; ++j ) {
    knots[deg+j] = a + h*j;
  }
}

double splines::omega(double x, size_t i, size_t ord, const std::vector<double> & knots) {
  double a = knots[i]; double b = knots[i+ord-1];
  if ( a < b ) {
    return (x-a)/(b-a);
  } else {
    return 0.0;
  }
}

double splines::basis(double x, size_t i, size_t ord, const std::vector<double> & knots) {
  if ( ord == 1 ) {
    double a = knots[i]; double b = knots[i+1];
    if ( a <= x && x < b ) {
      return 1.0;
    } else {
      return 0.0;
    }
  } else {
    double om0 = omega(x, i, ord, knots);
    double om1 = omega(x, i+1, ord, knots);
    double B0 = basis(x, i, ord-1, knots);
    double B1 = basis(x, i+1, ord-1, knots);
    return om0 * B0 + (1-om1) * B1;
  }
}

double splines::spline(double x, size_t ord, const std::vector<double> & knots,
      const std::vector<double> & weights) {
  double y = 0.0;
  for ( size_t i = 0; i < weights.size(); ++i ) {
    y += basis(x, i, ord, knots) * weights[i];
  }
  return y;
}

double Spline::operator()(double x) const {
  return spline(x, order, knots, weights);
}

Spline splines::import_spline(std::string filename) {
  Spline spl;
  size_t num_break_points;
  std::ifstream file(filename.c_str());
  // read order
  if ( file.good() ) {
    std::string order_str;
    std::getline(file, order_str);
    std::stringstream(order_str) >> spl.order;
  } else {
    throw std::runtime_error("unable to read order from spline file");
  }
  // function to read a vector
  auto read_vec = [](const std::string & line) -> std::vector<double> {
    std::vector<double> vec;
    std::stringstream ss(line);
    std::string word;
    while ( ss.good() ) {
      std::getline(ss, word, ' ');
      vec.push_back(atof(word.c_str()));
    }
    return vec;
  };
  // read break points
  if ( file.good() ) {
    std::string bp_line;
    std::getline(file, bp_line);
    std::vector<double> break_points = read_vec(bp_line);
    num_break_points = break_points.size();
    gen_knots_ref(spl.order, break_points, spl.knots);
  } else {
    throw std::runtime_error("unable to read break points from spline file");
  }
  // read weights
  if ( file.good() ) {
    std::string w_line;
    std::getline(file, w_line);
    spl.weights = read_vec(w_line);
  } else {
    throw std::runtime_error("unable to read weights from spline file");
  }
  // check that the number of weights is correct
  if ( spl.weights.size() != num_basis(spl.order, num_break_points) ) {
    throw std::runtime_error("incompatible order, weihts and break points");
  }
  return spl;
}


/// methods for piecewise linear function

inline double linfun(double x, double x0, double x1, double y0, double y1) {
  return y0 + (y1-y0)/(x1-x0) * (x-x0);
}

double PiecewiseLinFunc::operator()(double x) const {
  size_t n = break_points.size();
  if ( n < 2 || n != values_break_points.size() ) {
    std::logic_error("invalid definition of piecewise linear function" + RIGHT_HERE);
  }
  double x0, x1, y0, y1;
  if ( x < break_points.front() ) { // extrapolate
    y0 = values_break_points[0];
    y1 = values_break_points[1];
    x0 = break_points[0];
    x1 = break_points[1];
    return linfun(x, x0, x1, y0, y1);
  } // else interpolate
  // find interval
  for ( size_t i = 1; i < n; ++i ) {
    if ( x < break_points[i] ) {
      y0 = values_break_points[i-1];
      y1 = values_break_points[i];
      x0 = break_points[i-1];
      x1 = break_points[i];
      return linfun(x, x0, x1, y0, y1);
    }
  } // else, x >= all breakpoints: extrapolate!
  y0 = values_break_points[n-2];
  y1 = values_break_points[n-1];
  x0 = break_points[n-2];
  x1 = break_points[n-1];
  return linfun(x, x0, x1, y0, y1);
}


PiecewiseLinFunc splines::import_piecewise_lin_fun(std::string filename) {
  std::ifstream file(filename.c_str());
  // call more general function for istream input
  return import_piecewise_lin_fun(file);
}

std::list<PiecewiseLinFunc> splines::import_multiple_piecewise_lin_funs(std::string filename) {
  std::ifstream file(filename.c_str());
  std::list<PiecewiseLinFunc> funlist;
  while ( file.good() ) {
    try {
      PiecewiseLinFunc fun = import_piecewise_lin_fun(file);
      funlist.push_back(fun);
    } catch ( std::runtime_error & e ) {
      // ignore error
    }
  }
  return funlist;
}


PiecewiseLinFunc splines::import_piecewise_lin_fun(std::istream & in) {
  PiecewiseLinFunc fun;
  // function to read a vector
  auto read_vec = [](const std::string & line) -> std::vector<double> {
    std::vector<double> vec;
    std::stringstream ss(line);
    std::string word;
    while ( ss.good() ) {
      std::getline(ss, word, ' ');
      vec.push_back(atof(word.c_str()));
    }
    return vec;
  };
  // read break points
  if ( in.good() ) {
    std::string bp_line;
    std::getline(in, bp_line);
    fun.break_points = read_vec(bp_line);
  } else {
    throw std::runtime_error("unable to read break points from spline file");
  }
  // read weights
  if ( in.good() ) {
    std::string val_line;
    std::getline(in, val_line);
    fun.values_break_points = read_vec(val_line);
  } else {
    throw std::runtime_error("unable to read weights from spline file");
  }
  // check that the number of weights is correct
  if ( fun.values_break_points.size() != fun.break_points.size() ) {
    throw std::runtime_error("incompatible order, weihts and break points");
  }
  return fun;
}
