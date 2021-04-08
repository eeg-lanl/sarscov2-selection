#include "linalg.hpp"

#include <cmath>
#include <iostream>
#include <algorithm> // any_of
#include <numeric> // reduce
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>

#include "macros.hpp"

std::vector<double> calc_initial_condition(double xi, double beta, double gamma,
    double alpha, double nu, double omega) {
  // create jacobian matrix
  size_t n = 3;
  std::vector y0(3, 0.0);
	gsl_matrix* Jac = gsl_matrix_calloc(n, n); // MALLOC

  /* Set elements of Jac
   *
   *       [-alpha   beta*(1-xi)   0    ]
   * Jac = [ alpha  -(gamma+nu)    0    ]
   *       [ 0       nu           -omega]
   *
   * The zero elements are already set.
   */
  gsl_matrix_set(Jac, 0, 0, -alpha);
  gsl_matrix_set(Jac, 0, 1, beta*(1-xi));
  gsl_matrix_set(Jac, 1, 0, alpha);
  gsl_matrix_set(Jac, 1, 1, -(gamma+nu));
  gsl_matrix_set(Jac, 2, 1, nu);
  gsl_matrix_set(Jac, 2, 2, -omega);

	// compute eigenspace and select right eigenvector
	gsl_vector_complex* eval = gsl_vector_complex_alloc(n); // MALLOC
	gsl_matrix_complex* evec = gsl_matrix_complex_alloc(n, n); // MALLOC
	gsl_eigen_nonsymmv_workspace* w = gsl_eigen_nonsymmv_alloc(n); // MALLOC
	gsl_eigen_nonsymmv(Jac, eval, evec, w);
	gsl_vector_view re_eval = gsl_vector_complex_real(eval);
	size_t idx = gsl_vector_max_index(&re_eval.vector);
	//double r0 = gsl_vector_get(&re_eval.vector, idx);
	gsl_vector_complex_view h_complex = gsl_matrix_complex_column(evec, idx);
	gsl_vector_view h = gsl_vector_complex_real(&h_complex.vector);

	// store the initial state in vector y0
	gsl_vector_view y0_vec = gsl_vector_view_array(y0.data(), n);
	gsl_vector_memcpy(&y0_vec.vector, &h.vector);

  // free jacobian matrix and aux matrices, vectors and memory
  gsl_matrix_free(Jac); // FREE
  gsl_eigen_nonsymmv_free(w); // FREE
  gsl_vector_complex_free(eval); // FREE
  gsl_matrix_complex_free(evec); // FREE

	// re-scale and test
  double sum = std::reduce(y0.begin(), y0.end(), 0.0);

  if ( sum == 0.0 ) {
    std::cerr << "WARNING! sum of eigenvector is zero" + RIGHT_HERE << std::endl;
    // replace with uniform vector
    y0 = {0.5, 0.5, 0.0}; // fall-back vector for when all else fails.
    return y0;
  } // else, divide by sum

  std::for_each(y0.begin(), y0.end(), [sum](double & x){x /= sum;});

  // check that all elements are positive
  if ( std::any_of(y0.begin(), y0.end(), [](double x){return x < 0;}) ) {
    std::cerr << "WARNING! unable to make non-negative eigenvector" + RIGHT_HERE << std::endl;
    // replace with uniform vector
    y0 = {0.5, 0.5, 0.0}; // fall-back vector for when all else fails.
    return y0;
  }

  // else, all is well!
	return y0;
}
