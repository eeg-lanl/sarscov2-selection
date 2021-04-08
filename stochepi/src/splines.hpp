#ifndef SPLINES_HPP_
#define SPLINES_HPP_

#include <iostream>
#include <vector>
#include <list>

namespace splines {
  size_t num_knots(size_t order, size_t num_break_points);
  size_t num_basis(size_t order, size_t num_break_points);
  void gen_knots_ref(size_t order, const std::vector<double> & break_points, std::vector<double> & knots);
  void gen_unif_knots_ref(size_t order, size_t num_break_points, double a, double b, std::vector<double> & knots);
  double omega(double x, size_t i, size_t ord, const std::vector<double> & knots);
  double basis(double x, size_t i, size_t ord, const std::vector<double> & knots);
  double spline(double x, size_t ord, const std::vector<double> & knots, const std::vector<double> & weights);

  struct Spline {
    std::vector<double> knots;
    std::vector<double> weights;
    size_t order;
    double operator()(double x) const;
  };

  /** a spline file consists of 3 lines.
   * the first line gives the order
   * the second line gives the break points
   * the third line gives the weights
   * the break points and weights are separated by spaces
   */
  Spline import_spline(std::string filename);

  /** sometimes it is simpler to use a continuous, piecewise linear function
   * instead of a spline */
  struct PiecewiseLinFunc {
    std::vector<double> break_points;
    std::vector<double> values_break_points;
    double operator()(double x) const;
  };

  /** a linear function file is given by two
   * rows. The first row contains the breakpoints,
   * the second row contains the values at these break points.
   * (This method calls the stream method below)
   */
  PiecewiseLinFunc import_piecewise_lin_fun(std::string filename);
  /** use a stream instead of a filename */
  PiecewiseLinFunc import_piecewise_lin_fun(std::istream & in);
  /** import multiple pievewise linear functions concatenated in one file */
  std::list<PiecewiseLinFunc> import_multiple_piecewise_lin_funs(std::string filename);
}

#endif
