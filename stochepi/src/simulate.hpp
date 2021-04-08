#ifndef SIMULATE_HPP_
#define SIMULATE_HPP_

/** @file simulate.hpp
 * @short Ensamble simulations of a Model
 */

#include "model.hpp"

/** simulate function that opens a file with the given name
 * and then calls the more general simulate function
 */
void simulate(const Model & model, double tmax, double dt, size_t J,
      std::string filename, int threads, unsigned long seed);

/** simulate function that writes output to a general output stream object
 */
void simulate(const Model & model, double tmax, double dt, size_t J,
      std::ostream & stream, int threads, unsigned long seed);


#endif
