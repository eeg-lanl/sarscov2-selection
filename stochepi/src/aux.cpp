#include "aux.hpp"

#include <algorithm> // max_element
#include <stdexcept>
#include <list> // easier to sort...
#include <fstream> // reading files
#include <sstream> // parse strings

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>

#include "macros.hpp"

/* functions defined on interfaces */

/* Method for Indexable class */

long Indexable::getIndex() const {
	return index;
}
void Indexable::assignIndex(long index) {
	this->index = index;
}




/* methods for the RandomIndexable class */

int RandomIndexable::getRandomIndex() const { return random_index; }
void RandomIndexable::setRandomIndex(Rng & rng) { random_index = rng.integer(); }

bool compareByRandomIndex(RandomIndexable* left, RandomIndexable* right) {
	if ( left != nullptr && right != nullptr ) {
		return left->getRandomIndex() < right->getRandomIndex();
	}	else {
		if ( right == nullptr ) return true;
		else return false;
	}
}

/* methods for the Printable class */

std::ostream & operator<<(std::ostream & os, const Printable & obj) {
	obj.print(os);
	return os;
}


// auxiliary functions...

/* NB: inline functions and templates are defined in the header.
 */

double weighted_geometric_average(double a, double b, double w) {
	if ( a < 0.0 || b < 0.0 || w < 0.0 || w > 1.0 ) {
		throw std::logic_error("ERROR: invalid parameters" + RIGHT_HERE);
	}
	return pow(a, w) * pow(b, 1.0-w);
}

std::pair<Vec, bool> eigen_values(double a, double b, double c, double d) {
	double det = a*d - b*c;
	double tr = a+d;
	double D = tr*tr - 4*det;
	if ( D >= 0 ) {
		double sqrtD = sqrt(D);
		auto lambda = std::make_pair(0.5*(tr+sqrtD), 0.5*(tr-sqrtD));
		return std::make_pair(lambda, true); // true mean real
	} else {
		double sqrtAbsD = sqrt(abs(D));
		auto reimag = std::make_pair(0.5*tr, 0.5*sqrtAbsD);
		return std::make_pair(reimag, false); // false means complex (and not real)
	}
}

std::pair<Vec, bool> eigen_values(const Mat & M) {
	return eigen_values(M.a, M.b, M.c, M.d);
}

Vec eigen_vector_norm1(double a, double b, double c, double d, double lambda) {
	WARN_UNTESTED_FUN
	double denom = sqrt((lambda-a)*(lambda-a) + b*b);
	double u1 = b / denom;
	double u2 = (lambda - a) / denom;
	return std::make_pair(u1, u2);
}

Vec eigen_vector_norm1(const Mat & M, double lambda) {
	return eigen_vector_norm1(M.a, M.b, M.c, M.d, lambda);
}

Vec eigen_vector_sum1(double a, double b, double c, double d, double lambda) {
	double u1 = -b / (a - lambda - b);
	double u2 = (a - lambda) / (a - lambda - b);
	return std::make_pair(u1, u2);
}

Vec eigen_vector_sum1(const Mat & M, double lambda) {
	return eigen_vector_sum1(M.a, M.b, M.c, M.d, lambda);
}


// numerically stable log_sum_exp



double log_sum_exp(const std::vector<double> & xs) {
	/* we use the fact that
	 * log_sum_exp(x_1, x_2, x_3, ... x_n) =
	 *    log_sum_exp(x_n, log_sum_exp(x_{n-1}, .., x_1)...)
	 * and sort the list (x_i)_i for better numerical precision
	 */
	if ( xs.empty() ) {
		throw std::invalid_argument("vector cannot be empty" + RIGHT_HERE);
	} // else ...
	// copy the vector to sort it
	std::list<double> xs_sorted(xs.begin(), xs.end());
	xs_sorted.sort();
	double lse = xs_sorted.front();
	for ( auto it = ++xs_sorted.begin(); it != xs_sorted.end(); ++it ) {
		lse = log_sum_exp(lse, *it);
	}
	return lse;
}



// inverse simpson index serves as an effective sample size
double log_eff_sam_size(const std::vector<double> & loglikes) {
	// log(2D) = log(sum(exp(l))^2 / sum(exp(2l)))
  if ( loglikes.empty() ) {
		throw std::invalid_argument("vector cannot be empty" + RIGHT_HERE);
  } // else...
  std::list<double> lls(loglikes.begin(), loglikes.end());
  lls.sort();
  double lse = lls.front();
  double lse2 = 2*lls.front();
  for ( auto it = ++lls.begin(); it != lls.end(); ++it ) {
    double & ll = *it;
    lse = log_sum_exp(lse, ll);
    lse2 = log_sum_exp(lse2, 2*ll);
  }
  return 2*lse - lse2;
}




/** Simple function to import a matrix from a csv file */
std::vector<std::vector<double>> import_mat(const std::string & filename) {
  std::vector<std::vector<double>> mat;
  std::ifstream file(filename.c_str());
  while ( file.good() ) {
    std::string line;
    std::vector<double> row;
    std::getline(file, line);
		if ( line.empty() ) {
			continue;
		}
    std::stringstream ss(line);
    std::string word;
    while ( ss.good() ) {
      std::getline(ss, word, ',');
      row.push_back(atof(word.c_str()));
    }
    mat.push_back(row);
  }
  return mat;
}


/** Simple function to import a list of numbers from a file.
 * Elements of the list must be separated by newlines ('\n')
 */
std::vector<double> import_vec(const std::string & filename) {
  std::vector<double> list;
  std::ifstream file(filename.c_str());
  while ( file.good() ) {
    std::string line;
    std::getline(file, line);
		if ( line.empty() ) {
			continue;
		}
    list.push_back(atof(line.c_str()));
  }
  return list;
}



std::pair<double, bool> spectral_radius(const std::vector<std::vector<double>> & mat) {
	size_t n = mat.size();
	if ( n == 0 ) {
		return std::make_pair(0.0, false);
	}

	// allocate space...
	gsl_eigen_nonsymm_workspace* workspace = gsl_eigen_nonsymm_alloc(n);
	gsl_matrix* A = gsl_matrix_alloc(n, n);
	gsl_vector_complex* eval = gsl_vector_complex_alloc(n);

	// copy mat to A
	for ( size_t i = 0; i < n; ++i ) {
		if ( mat[i].size() != n ) {
			// matrix is not square
			return std::make_pair(0.0, false);
		}
		for ( size_t j = 0; j < n; ++ j) {
			gsl_matrix_set(A, i, j, mat[i][j]);
		}
	}

	// compute eigenvalues
	int ok = gsl_eigen_nonsymm(A, eval, workspace);

	if ( ok != GSL_SUCCESS ) {
		return std::make_pair(0.0, false);
	}

	// find the largest eigenvalue
	double max_eval = 0.0;
	bool init_max_eval = false;
	for ( size_t i = 0; i < n; ++i ) {
		double x = gsl_complex_abs(gsl_vector_complex_get(eval, i));
		if ( !init_max_eval || max_eval < x ) {
			max_eval = x; init_max_eval = true;
		}
	}

	// C-style cleaunp required for GSL routines
	gsl_vector_complex_free(eval);
	gsl_matrix_free(A);
	gsl_eigen_nonsymm_free(workspace);

	// return spectral radius
	return std::make_pair(max_eval, true);
}
