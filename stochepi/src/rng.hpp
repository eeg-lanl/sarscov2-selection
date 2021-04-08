#ifndef RNGCLASS_HPP
#define RNGCLASS_HPP

#include <vector>

#include <gsl/gsl_rng.h>


/** this is just a c++ class wrapper for GNU's gsl_rng and it's
 * gsl_randist, except for (at the moment) the methods Erlang, PosBinomial,
 * SkewNormal, ...
 * Also, some parameter conventions might be different
 * such as in BetaMnSd and Exponential.
 * TODO: replace with C++11 alternative? At least, make "compatible" with
 * C++11 <random> header.
 */
class Rng {
public:
	Rng();
	Rng(unsigned long );
	Rng(const Rng & ) = delete; // copy constructor not defined
	void operator=(const Rng & ) = delete; // assignment operator not defined
	void seed(unsigned long ); // seed with 32 bit integer
	virtual ~Rng(); // destructor
// members to make Rng compatible with the <random> library (and e.g. std::shuffle)
	typedef unsigned long result_type;
  result_type operator()();
  result_type min() const;
  result_type max() const;
// distributions
	unsigned long integer(); // call rand_int32
	unsigned long integer(unsigned long ); // modulo max
	bool bit(); // a random bit
	bool bernoulli(double ); // Bernoulli(p) with values 0 = false, 1 = true
	double uniform(); // uniform(0,1)
	double uniform(double , double ); // Uniform(a,b)
	double normal(double , double ); // mean and sd
	double stdNormal(); // normal(0, 1)
	double logNormal(double , double ); // location, scale parameterization
	double logNormalMnSd(double , double ); // Mean, Sd parameterization
	double skewNormal(double , double , double ); // mu, sigma and alpha
	double exponential(double ); // pass the rate parameter (1/mean)
	double truncExponential(double , double); // less than second argument, naive rejection method (todo)
	double weibull(double , double ); // scale, shape
	double erlang(double , unsigned ); // simple Erlang(rate,shape) (sum of exponentials)
	double gamma(double , double ); // Gamma(scale,shape) distribution
	double paretoMinusOne(double ); // Pareto(a, 1)-1
	double beta(double , double ); // Beta(alpha,beta) distribution
	double betaMnSd(double , double ); // Beta, but with mean and sd as parameterization
	int binomial(int , double ); // Binom(n,p)
	int posBinomial(int , double); // Binomial, conditioned on > 0
	int poisson(double ); // Poisson(lambda) with mean lambda
 	int posPoisson(double ); // a.k.a. Zero-truncated Poisson distribution
	int zeta(double ); // the Zipf or zeta distribution (power law)
	int hypergeometric(int , int , int ); // Hyper(n, K, N) (N = population size, K = successes, n = draws)
  int skellam(double , double ); // Skellam(lambda_1, lambda_2) = Poisson(lambda_1) - Poisson(lambda_2)
  int symSkellam(int, double ); // Symmetric Skellam distribution with location parameter and sd
	void shuffle(int* , int , int ); // a random sample (unique)
  void dirichlet(double* , int , const double* ); // the Dirichlet distribution
	int categorical(const std::vector<double> & ); // vector of probabilities
	/* TODO cathegorical becomes MUCH more efficient when the size of the sample is large
	 * Implement a version with multiple samples
	 */
private:
	gsl_rng* gslRng;
	void init();
};


#endif /* RNGCLASS_HPP_ */
