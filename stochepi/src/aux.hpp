#ifndef AUX_HPP_
#define AUX_HPP_

#include <cmath> // log, exp
#include <string>
#include <vector>
#include <complex>
#include <gsl/gsl_math.h> // log_1px

#include "rng.hpp"
#include "macros.hpp"

/* typedefs */

typedef double Real;
typedef std::complex<Real> Complex;
typedef std::vector<Real> RealVec;
typedef std::vector<Complex> ComplexVec;

/* enums */

enum CensorCode {
  UNCENSORED_CODE = 0,
  LEFT_CENSORED_CODE = 1,
  RIGHT_CENSORED_CODE = 2,
  MISSING_CODE = 3,
  INTERVAL_CENSORED_CODE = 4
};

/* constants */

constexpr double UL_PER_ML = 1e3;

/* interfaces */

class Printable {
public:
	Printable() { /* empty */ }
	virtual ~Printable() { /* empty */ }
	virtual void print(std::ostream & ) const = 0;
};

std::ostream & operator<<(std::ostream & , const Printable & );

/** a class for giving list idems a index that can be different from
 * the order of the list. Defining invalid indices (e.g. negative integers)
 * allows for skipping of items
 */
class Indexable {
public:
	const static long invalid_index = -1;
	Indexable() : index(invalid_index) {}
	virtual ~Indexable() { /* empty */ }
  long getIndex() const;
  void assignIndex(long );
protected:
  long index;
};


/** a class that allows for random shuffling by sorting a list
 * For instant std::list<Agent*> or std::list<Contact*>
 */
class RandomIndexable {
public:
	RandomIndexable() : random_index(0) {}
	virtual ~RandomIndexable() {}
	int getRandomIndex() const;
	// use an RNG to reset the random_index
	void setRandomIndex(Rng & );
private:
 	unsigned long random_index;
};

bool compareByRandomIndex(RandomIndexable* , RandomIndexable* );


/** a class that can be inherited instead of Base
 * if the dup method needs to be inherited from
 * the base class.
 *
 * @todo is there a way to mimick covariant return types for the dup method?
 *
 * Example usage:
 * ```
 * class Base { public virtual Base* dup() const; ... };
 * // Base must declare dup (could be pure)
 * class Derived : public Cloneable<Base, Derived> { ... };
 * // Derived now has a dup method that clones the correct object
 * ```
 */
template<class Base, class Derived>
class Cloneable : public Base {
public:
  virtual ~Cloneable();
  virtual Base* dup() const override { // TODO: can this returen Derived???
    return new Derived(*static_cast<const Derived*>(this));
  }
};



/* inline functions */

inline double logit(double x) { return log(x/(1-x)); }
inline double expit(double x) { // inverse of logit
  if ( x > 0 ) { // numerically stable method..
    return 1/(1+exp(-x));
  } else {
    double y = exp(x);
  	return y/(1+y);
  }
}

inline double log_sum_exp(double x, double y) {
	/* log(exp(x) + exp(y)) = x + log(1 + exp(y-x)) */
	if ( x > y ) {
		return x + gsl_log1p(exp(y-x));
	} else {
		return y + gsl_log1p(exp(x-y));
	}
}

inline double log_diff_exp(double x, double y) {
	if ( x <= y ) {
		throw std::invalid_argument("x must be larger than y" + RIGHT_HERE);
	}
	// x > y, hence y-x < 0, hence exp(y-x) < 1
	return x + gsl_log1p(-exp(y-x));
}

inline double extinction_threshold(double u, double x0, double R0) {
    /* calculate the extinction threshold for a declining population (R0 < 1).
     * The threshold for switching between ODEs and a stochastic description
     * is at x0, and u must be a Uniform(0,1) deviate
     */
    double squ = pow(u, 1.0/x0);
    return x0*(1-squ)/(1-squ*R0);
}
inline double major_outbreak_prob(double R) {
    /* The probability of a "major outbreak" given a reproduction number.
     * this is for adding mutants and activating effectors.
     * Notice that this function has a cusp at R = 1
     */
    return ( R > 1.0 ? 1-1/R : 0.0 );
}
inline double deriv_major_outbreak_prob(double R) {
    /* the derivative of major_ourbreak_prob
     * Notice that 1-1/R is not differentiable at R = 1
     */
    return ( R > 1.0 ? 1/(R*R) : 0.0 );
}
inline double major_outbreak_prob_C1(double R, double epsilon=0.1) {
    /* a C^1 approximation of major_outbreak_prob given by
     * (R-1)^2/(R*(R-1) + epsilon)
     * Notice that when epsilon = 0, we get (R-1)^2/(R*(R-1)) = 1 - 1/R
     */
    return ( R > 1.0 ? (R-1)*(R-1) / (R*(R-1) + epsilon) : 0.0 );
}
inline double deriv_major_outbreak_prob_C1(double R, double epsilon=0.1) {
    // the derivative (w.r.t. R) of major_outbreak_prob_C1
    double nmr = R*(R-1) + epsilon;
    return (R-1)*(2*epsilon+R-1)/(nmr*nmr);
}
inline int mod(int m, int n) {
    // compute m mod n with a representative 0 <= r <= n
    if ( m < 0 ) m += ((-m)/n + 1)*n;
    return m % n;
}

template<class T>
T quick_power(T x, int n) {
	// computes x^n by repeatedly squaring x. by definition, 0^0 = 1
	if ( n < 0 ) {
		x = 1.0/x;
		n = -n;
	}
  T result = 1;
  while ( n ) {
    if ( n & 1 ) {
			result *= x;
		}
    n >>= 1;
    x *= x;
  }
  return result;
}

template<class T>
T monomial(const std::vector<T> & x, const std::vector<int> & n) {
  T z = 1.0;
  for ( size_t i = 0; i < std::min(n.size(), x.size()); ++i ) {
    z *= quick_power(x[i], n[i]);
  }
  return z;
}


/* return a weighted geometric average of a and b,
 * where w is the weight of a
 */
double weighted_geometric_average(double a, double b, double w);


// alias for std::pair<double, double>
typedef std::pair<double, double> Vec;

// a 2x2 matrix
struct Mat {
	Mat() : a(0), b(0), c(0), d(0) {}
	Mat(double a, double b, double c, double d) : a(a), b(b), c(c), d(d) {}
	double a, b, c, d;
};

/* compute the eigen values of a 2x2 matrix. (a, b; c, d)
 * second argument of the pair is true if the eigenvalues are real,
 * and false if the eigenvalues are complex.
 * if the eigenvalues are real, the first element of the pair
 * consists of the two eigenvalues.
 * if the eigenvalues are complex, the real and imaginary part are returned
 */
std::pair<Vec, bool> eigen_values(double a, double b, double c, double d);
std::pair<Vec, bool> eigen_values(const Mat & M);

/* compute the eigen vector of (a, b; c, d) corresponding to eigenvalue lambda
 * this function assumes that the eigenvalue is correct
 * eigen_vector_norm1 is a vector of Euclidean norm.
 * eigen_vector_sum1 is a vector that sums to 1
 */
Vec eigen_vector_norm1(double a, double b, double c, double d, double lambda);
Vec eigen_vector_norm1(const Mat & M, double lambda);

Vec	eigen_vector_sum1(double a, double b, double c, double d, double lambda);
Vec	eigen_vector_sum1(const Mat & M, double lambda);

// log_sum_exp for lists and vectors


double log_sum_exp(const std::vector<double> & xs);

double log_eff_sam_size(const std::vector<double> & loglikes);


std::vector<std::vector<double>> import_mat(const std::string & filename);
std::vector<double> import_vec(const std::string & filename);


std::pair<double, bool> spectral_radius(const std::vector<std::vector<double>> & mat);


#endif
