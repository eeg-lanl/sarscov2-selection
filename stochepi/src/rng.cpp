#include "rng.hpp"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_zeta.h>
#include <stdexcept>
#include <cmath>
#include <climits>
#include <numeric> // partial sum
#include <sstream> // formatting error messages

#include "macros.hpp"



/******* members for Rng ********/

Rng::Rng() {
	init(); // uses standard seed
}
Rng::Rng(unsigned long s) {
	init();	seed(s);
}
Rng::~Rng() {
	gsl_rng_free(gslRng);
}
void Rng::init() {
	gslRng = gsl_rng_alloc(gsl_rng_taus);
}
void Rng::seed(unsigned long s) {
	gsl_rng_set(gslRng, s);
}

// compatibility with <random>
Rng::result_type Rng::operator()() {
    return integer();
}
Rng::result_type Rng::min() const {
	return 0;
}
Rng::result_type Rng::max() const {
	return gsl_rng_max(gslRng) - gsl_rng_min(gslRng);
}
// distributions

unsigned long Rng::integer() {
	return gsl_rng_get(gslRng) - gsl_rng_min(gslRng);
}
unsigned long Rng::integer(unsigned long max) { // modulo max
	return gsl_rng_uniform_int(gslRng, max);
}
bool Rng::bit() {
	return gsl_rng_uniform_int(gslRng, 2) == 1;
}
bool Rng::bernoulli(double p) {
	if ( p < 0.0 || p > 1.0 ) {
		throw std::invalid_argument("p must be between 0 and 1" + RIGHT_HERE);
	}
	return gsl_ran_bernoulli(gslRng, p);
}

double Rng::uniform() { // uniform(0,1), excludes 0 and 1
	 return gsl_ran_flat(gslRng, 0.0, 1.0);
}

double Rng::uniform(double a, double b) { // Uniform(a,b)
	 if ( a <= b ) {
		 if ( a==b ) {
			 return a;
		 } else {
			 return gsl_ran_flat(gslRng, a, b);
		 }
	 } else {
		 std::stringstream ss;
		 ss << "a (" << a << ") must be smaller than b (" << b << ")" << RIGHT_HERE;
		 throw std::invalid_argument(ss.str());
	 }
}

double Rng::normal(double mu, double sigma) { // mean and sd
	if ( sigma >= 0.0 ) {
		if ( sigma == 0.0 ) {
			return mu;
		}	else {
			return mu + gsl_ran_gaussian(gslRng, sigma);
		}
	}	else {
		throw std::invalid_argument("sigma must be non-negative" + RIGHT_HERE);
	}
}

double Rng::stdNormal() {
	return gsl_ran_gaussian(gslRng, 1.0);
}

double Rng::logNormal(double m, double s) { // location, scale parameterization
	if ( s < 0.0 ) {
		throw std::invalid_argument("s must be non-negative" + RIGHT_HERE);
	}
	if ( s == 0.0 ) return exp(m);
	return exp(normal(m, s));
}

double Rng::logNormalMnSd(double mu, double sigma) { // Mean, Variance parameterization
	if ( mu <= 0.0 || sigma < 0.0 ) {
		throw std::invalid_argument("non-positive mu or negative sigma" + RIGHT_HERE);
	}
	if ( sigma == 0.0 ) return mu;
	// now handle non-trivial cases...
	double x = 1.0 + (sigma*sigma) / (mu*mu);
	double s = sqrt(log(x));
	double m = log(mu/sqrt(x));
	return logNormal(m, s);
}

double Rng::skewNormal(double m, double s, double alpha) {
	if ( s < 0.0 ) {
		throw std::invalid_argument("s must be nonnegative" + RIGHT_HERE);
	}
	if ( s == 0.0 ) return m; // in the limit s --> 0
	// now handle non-trivial cases
	double Z1 = normal(0.0,1.0);
	if ( alpha == 0.0 ) return s * Z1 + m; // not skewed, just normal
	else {
		// X = 1/sqrt{1+alpha^2} alpha |Z1| + Z2
		double Z2 = normal(0.0,1.0);
		double absZ1 = (Z1 > 0 ? Z1 : -Z1);
		return m + s * ((alpha*absZ1 + Z2) / sqrt(1.0 + alpha*alpha));
	}
}

double Rng::exponential(double lambda) {
	if ( lambda > 0 ) {
		return gsl_ran_exponential(gslRng, 1.0/lambda);
	} else {
		throw std::invalid_argument("lambda must be positive" + RIGHT_HERE);
	}
}
double Rng::truncExponential(double lambda, double max) {
	double x;
	do {
		x = exponential(lambda);
	} while ( x > max );
	return x;
}

double Rng::weibull(double a, double b) { // scale, shape
	if ( a >= 0.0 && b >= 0.0 ) {
		if ( a == 0 ) {
			return 0.0;
		}	else {
			return gsl_ran_weibull(gslRng, a, b);
		}
	}	else {
		throw std::invalid_argument("a and b must be non-negative" + RIGHT_HERE);
	}
}

double Rng::erlang(double lambda, unsigned k) { // simple Erlang(rate,shape) (sum of exponentials)
	if ( lambda > 0 ) {
		double ans = 0.0;
		for ( unsigned i = 0; i < k; ++i ) {
			ans += exponential(lambda);
		}
		return ans;
	} else {
		throw std::invalid_argument("lambda must be positive" + RIGHT_HERE);
	}
}

double Rng::gamma(double scale, double shape) { // Gamma(scale,shape) distribution {
	if ( scale < 0.0 || shape <= 0.0 ) {
		throw std::invalid_argument("negative scale or nonnegative shape" + RIGHT_HERE);
	}
	if ( scale == 0.0 ) {
		return 0.0;
	}	else {
		return gsl_ran_gamma(gslRng, shape, scale);
	}
}

double Rng::paretoMinusOne(double a) {
	if ( a <= 0 ) {
		throw std::invalid_argument("a must be positive" + RIGHT_HERE);
	}
	return gsl_ran_pareto(gslRng, a, 1.0) - 1.0;
}

double Rng::beta(double alpha, double beta) { // Beta(alpha,beta) distribution
	if ( alpha > 0.0 && beta > 0.0 ) {
		return gsl_ran_beta(gslRng, alpha, beta);
	}	else {
		throw std::invalid_argument("alpha and beta must be positive" + RIGHT_HERE);
	}
}

double Rng::betaMnSd(double M, double SD) {
	if ( M <= 0.0 || M >= 1.0 || SD <= 0.0 ) {
		throw std::invalid_argument("M must be between 0 and 1 and SD must be positive" + RIGHT_HERE);
	}
	double a, sumab;
	sumab = M*(1.0-M)/(SD*SD) - 1.0;
	a = M*sumab;
	return beta(a,sumab-a);
}

int Rng::binomial(int n, double p) { // Binom(n,p)
	if ( p >= 0.0 && p <= 1.0 && n >= 0 ) {
		return gsl_ran_binomial(gslRng, p, unsigned(n));
	}	else {
		throw std::invalid_argument("p must be between 0 and 1 and n must be non-negative" + RIGHT_HERE);
	}
}

int Rng::posBinomial(int n, double p) {
	/* gets a sample x from the binomial distribution
	 * conditioned on x > 0.
	 */
	// handle some trivial cases
	int x = 0;
	if ( n <= 0 || p <= 0.0 || p > 1.0 ) {
		throw std::invalid_argument("n must be positive and p must be between 0 and 1" + RIGHT_HERE);
	}
	if ( n == 1 ) return 1;
	// handle non-trivial cases
	double np = n*p;
	if ( np > 1 ) {
		// use a naive method for large np
		do {
			x = binomial(n,p);
		} while ( x == 0 );
	} else {
		if ( np > 1e-6 ) {
			// chop-down method
			int n1 = n; double qn_aux = 1.0-p; double qn = 1.0;
			while ( n1 ) { // fast method for calculating qn := (1-p)^n
				if (n1 & 1) qn *= qn_aux;
				qn_aux *= qn_aux;  n1 >>= 1;
			}
			double pp = p / (1-p);
			double pdf = n * pp * qn; // probability of x = 1
			int ell = 1; n1 = n;
			double ran = uniform() * (1.0-qn); // condition on x > 0
			while ( ran > pdf && ell < n ) {
				ran -= pdf;
				ell++; n1--;
				pdf *= pp * n1;
				ran *= ell;
			}
			x = ell;
		}	else {
			/* use Poisson method
			 * First make an exp(1) jump, and scale so that it is between 0 and np. Then start
			 * counting exp(1) jumps below np.
			 */
			double u0 = uniform(); double uell(0.0);
			double expr_np = exp(-np); // don't take logs all the time...
			double expr_j0 = 1.0-u0*(1.0-expr_np); // take products, calc exp(np) and don't calc logs!!
			if ( expr_j0 <= expr_np ) {
				x = 1;
			}	else {
				int ell = 0;
				double expr_jell = expr_j0;
				while ( expr_jell > expr_np ) {
					ell++; // take another exponential jump.
					uell = uniform();
					expr_jell *= uell;
				}
				x = ell;
			}
		}
	}
	return x;
}

int Rng::poisson(double lambda) { // Poisson(lambda)
	if ( lambda >= 0.0 ) {
		if ( lambda == 0.0 ) {
			return 0;
		}	else {
			return gsl_ran_poisson(gslRng, lambda);
		}
	}	else {
		throw std::invalid_argument("lambda must be non-negative" + RIGHT_HERE);
	}
}

int Rng::posPoisson(double lambda) {
  if ( lambda < 0.0 ) {
		throw std::invalid_argument("lambda must be non-negative" + RIGHT_HERE);
	}
  int x = 0.0; // return value
  if ( lambda == 0.0 ) {
    return x;
  } else {
    if ( true ) { // TODO: case for large lambda
    	// Poisson with rejection
    	do {
    		x = gsl_ran_poisson(gslRng, lambda);
    	} while ( x == 0 );
    	return x;
  	}
		// else TODO: case for small lambda
	}
}

int Rng::zeta(double s) {
	/* dumb implementation: sample a uniform(0,zeta(s)) number u
	 * and start taking n^{-s} until u < 0
	 * TODO: use pareto and a rejection scheme
	 */
	if ( s <= 1.0 ) {
		throw std::invalid_argument("s must be > 1" + RIGHT_HERE);
	}
	double u = uniform(0.0, gsl_sf_zeta(s));
	int n = 1;
	while ( true ) {
		u -= pow(n, -s);
		if ( u <= 0.0 ) break;
		else ++n;
	}
	return n;
}



int Rng::hypergeometric(int n, int K, int N) { // draws, successes, population size
	if ( N >= 0 && K >= 0 && n >= 0 && K <= N && n <= N ) {
		return gsl_ran_hypergeometric(gslRng, K, N-K, n);
	}	else {
		throw std::invalid_argument("invalid parameters" + RIGHT_HERE);
	}
}


int Rng::skellam(double lambda1, double lambda2) {
    if ( lambda1 < 0.0 || lambda2 < 0.0 ) {
			throw std::invalid_argument("lambda1 or lambda2 is negative" + RIGHT_HERE);
		}
    return poisson(lambda1) - poisson(lambda2);
}

int Rng::symSkellam(int loc, double sd) {
    if ( sd < 0 ) {
			throw std::invalid_argument("sd is negative" + RIGHT_HERE);
		}
    double lambda = 0.5*sd*sd;
    return loc + skellam(lambda, lambda);
}

void Rng::shuffle(int* xs, int min, int n) { // a random sample (unique)
	if ( n > 0 ) {
		for ( int i = 0; i < n; ++i ) {
			xs[i] = min + i;
		}
		gsl_ran_shuffle(gslRng, xs, size_t(n), sizeof(int));
	}	else {
		throw std::invalid_argument("n must be positive" + RIGHT_HERE);
	}
}

void Rng::dirichlet(double* xs, int n, const double* alphas) {
  if ( n > 0 ) {
    // TODO: accept zero alphas: sample is then restricted to subspace.
    gsl_ran_dirichlet(gslRng, n, alphas, xs);
  } else {
    throw std::invalid_argument("n must be positive" + RIGHT_HERE);
  }
}

int Rng::categorical(const std::vector<double> & p) { // number of categories and probabilities
	 gsl_ran_discrete_t* g = gsl_ran_discrete_preproc(p.size(), p.data());
	 int category = gsl_ran_discrete(gslRng, g);
	 gsl_ran_discrete_free(g);
	 return category;
}
