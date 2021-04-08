#include "distributions.hpp"

#include <cmath> // fabs
#include <exception>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h> // the log CDF and complementary CDF of Gaussian

#include "aux.hpp"
#include "macros.hpp"

/* simply use the poisson distribution for dispersion (d) smaller than
 * the following lower bound
 */
constexpr double POLYA_DISPERSION_LOWER_BOUND = 1e-4;



// likelihoods...

double ll_binomial(int n, double p, int k) {
	/* likelihood = {n choose k} * p^k * (1-p)^{n-k}
	 */
	return  gsl_sf_lnchoose(n, k) + k * log(p) + (n-k) * log(1-p);
}

/** likelihood = \binom{n, k} B(pr + k, (1-p)r + n-k) / B(pr, (1-p)r)
 */
double ll_beta_binomial(int n, double p, double r, int k) {
	double a = p*r;	double b = (1-p)*r;
	return gsl_sf_lnchoose(n, k) + gsl_sf_lnbeta(a + k, b + n - k) - gsl_sf_lnbeta(a, b);
}

double ll_poisson(double lambda, int k) {
	/* likelihood = lambda^k/k! * exp(-lambda),
	 * so log likelihood = k * log(lambda) - log(k!) - lambda
	 */
	return k * log(lambda) - lambda - gsl_sf_lngamma(k+1);
}

double ll_neg_binomial(double lambda, double n, int k) {
	/* likelihood is Gamma(n + k) / (Gamma(k+1) Gamma(n)) * p^n * (1-p)^k
	 * where p = lambda / (n + lambda) and n = 1/d
	 */
	double p = lambda / (n+lambda);
	return gsl_sf_lngamma(n+k) - gsl_sf_lngamma(k+1) - gsl_sf_lngamma(n) +
		n * log(1.0-p) + k * log(p);
}

double ll_polya(double lambda, double d, int k) {
	if ( d < POLYA_DISPERSION_LOWER_BOUND ) {
		return ll_poisson(lambda, k);
	}
	else {
		return ll_neg_binomial(lambda, 1/d, k);
	}
}


// normal PDF, CDF and CCDF

double norm_lpdf(double x, double mu, double sigma) {
  double z = (x-mu)/sigma;
  return -0.5 * z*z - log(M_SQRT2 * M_SQRTPI * sigma);
}

double norm_lcdf(double x, double mu, double sigma) {
  double z = (x-mu)/sigma;
  return gsl_sf_log_erfc(-z/M_SQRT2) - M_LN2;
}

double norm_lccdf(double x, double mu, double sigma) {
  double z = (x-mu)/sigma;
  return gsl_sf_log_erfc(z/M_SQRT2) - M_LN2;
}


/* Returns the log-likelihood, and the validity (always true)
 */
std::pair<double, bool> censored_norm_lpdf(double x, double xhat, double sigma, int c) {
	std::pair<double, bool> log_like(0.0, true); // to-be-returned
	switch ( c ) {
		case UNCENSORED_CODE: {
			log_like.first += norm_lpdf(x, xhat, sigma);
			break;
		}
		case LEFT_CENSORED_CODE: {
			log_like.first += norm_lcdf(x, xhat, sigma);
			break;
		}
		case RIGHT_CENSORED_CODE: {
			log_like.first += norm_lccdf(x, xhat, sigma);
			break;
		}
		case MISSING_CODE: {
			/* loglike += 0.0 */
			break;
		}
		case INTERVAL_CENSORED_CODE: {
			WARN_NOT_IMPLEMENTED
			break;
		}
		default: {
			std::stringstream ss;
			ss << "invalid censoring code '" << c << "' " <<  RIGHT_HERE;
			throw std::logic_error(ss.str());
			break;
		}
	}
	return log_like;
}




// a combination of the above (log-normal) that takes censoring into account

/** Returns the log-likelihood, and the validity.
 * i.e. if the likelihood is 0, the log-likelihood is -infinity,
 * and this is indicated by a 'false' return value
 * the data x MUST BE POSITIVE
 */
std::pair<double, bool> censored_lognorm_lpdf(double x, double xhat, double sigma, int c) {
	if ( x <= 0.0 && c != MISSING_CODE ) {
		throw std::logic_error("data x must be a positive value" + RIGHT_HERE);
	} // else...
  std::pair<double, bool> log_like(0.0, true); // to-be-returned
  switch ( c ) {
    case UNCENSORED_CODE: {
      if ( xhat > 0 ) {
        log_like.first += norm_lpdf(log(x), log(xhat), sigma);
      } else {
        log_like.second = false;
      }
      break;
    }
    case LEFT_CENSORED_CODE: {
      if ( xhat > 0 ) {
        log_like.first += norm_lcdf(log(x), log(xhat), sigma);
      } else {
        /* log_like += 0.0, because the detection limit
         * is infinitly far away from the actual value
         */
      }
      break;
    }
    case RIGHT_CENSORED_CODE: {
      if ( xhat > 0 ) {
        log_like.first += norm_lccdf(log(x), log(xhat), sigma);
      } else {
        log_like.second = false;
      }
      break;
    }
    case MISSING_CODE: {
      /* loglike += 0.0 */
      break;
    }
    case INTERVAL_CENSORED_CODE: {
      WARN_NOT_IMPLEMENTED
      break;
    }
    default: {
      std::stringstream ss;
      ss << "invalid censoring code '" << c << "' " <<  RIGHT_HERE;
      throw std::logic_error(ss.str());
      break;
    }
  }
  return log_like;
}



std::pair<double, bool> binomial_lpmf(int k, int n, double p) {
	// check validity of arguments
	if ( p < 0 || p > 1 ) {
		std::stringstream ss;
		ss << "probability p must be in [0,1], but is equal to " << p;
		throw std::invalid_argument(ss.str());
	}
	if ( n < 0 ) {
		std::stringstream ss;
		ss << "sample size n must be non-negative, but is equal to " << n;
		throw std::invalid_argument(ss.str());
	}
	// handle individual cases
	if ( p > 0 && p < 1 ) { // generic case
		if ( k < 0 || k > n ) {
			return std::make_pair(0.0, false);
		} else {
			double ll = ll_binomial(n, p, k);
			return std::make_pair(ll, true);
		}
	} else {
		if ( p == 0 ) {
			 return ( k == 0 ? std::make_pair(0.0, true) : std::make_pair(0.0, false) );
		} else { // p == 1
			return ( k == n ? std::make_pair(0.0, true) : std::make_pair(0.0, false) );
		}
	}
}

std::pair<double, bool> beta_binomial_lpmf(int k, int n, double p, double r) {
	// check validity of arguments
	if ( p < 0 || p > 1 ) {
		std::stringstream ss;
		ss << "probability p must be in [0,1], but is equal to " << p;
		throw std::invalid_argument(ss.str());
	}
	if ( n < 0 ) {
		std::stringstream ss;
		ss << "sample size n must be non-negative, but is equal to " << n;
		throw std::invalid_argument(ss.str());
	}
	if ( r <= 0 ) {
		std::stringstream ss;
		ss << "overdispersion parameter r must be positive, but is equal to " << r;
		throw std::invalid_argument(ss.str());
	}
	// handle individual cases
	if ( p > 0 && p < 1 ) { // generic case
		if ( k < 0 || k > n ) {
			return std::make_pair(0.0, false);
		} else {
			double ll = ll_beta_binomial(n, p, r, k);
			return std::make_pair(ll, true);
		}
	} else {
		if ( p == 0 ) {
			 return ( k == 0 ? std::make_pair(0.0, true) : std::make_pair(0.0, false) );
		} else { // p == 1
			return ( k == n ? std::make_pair(0.0, true) : std::make_pair(0.0, false) );
		}
	}
}


std::pair<double, bool> poisson_lpmf(int k, double mu) {
	if ( mu < 0 ) {
		std::stringstream ss;
		ss << "parameter mu must be non-negative, but is equal to " << mu;
		throw std::invalid_argument(ss.str());
	} // else...
	if ( mu > 0 ) { // generic case
		if ( k < 0 ) {
			return std::make_pair(0.0, false);
		} else {
			double ll = ll_poisson(mu, k);
			return std::make_pair(ll, true);
		}
	} else { // mu = 0: only possible mu k = 0
		return (k == 0 ? std::make_pair(0.0, true) : std::make_pair(0.0, false));
	}
}

std::pair<double, bool> neg_binomial_lpmf(int k, double mu, double r) {
	if ( mu < 0 ) {
		std::stringstream ss;
		ss << "parameter mu must be non-negative, but is equal to " << mu;
		throw std::invalid_argument(ss.str());
	} // else...
	if ( r <= 0 ) {
		std::stringstream ss;
		ss << "parameter r must be positive, but is equal to " << r;
		throw std::invalid_argument(ss.str());
	}
	if ( mu > 0 ) { // generic case
		if ( k < 0 ) {
			return std::make_pair(0.0, false);
		} else {
			double ll = ll_neg_binomial(mu, r, k);
			return std::make_pair(ll, true);
		}
	} else { // mu = 0: only possible mu k = 0
		return (k == 0 ? std::make_pair(0.0, true) : std::make_pair(0.0, false));
	}
}



/* methods for priors */

Prior::Prior() { /* empty */ }
Prior::~Prior() { /* empty */ }
Prior* Prior::dup() const {
	return new Prior(*this);
}
double Prior::loglike(double x) const {
	return 0.0;
}

UniformPrior::UniformPrior(double a, double b) : a(a), b(b) {
	if ( a >= b ) {
		throw std::invalid_argument("a must be smaller than b" + RIGHT_HERE);
	}
	ll = -log(b - a);
}
UniformPrior* UniformPrior::dup() const {
	return new UniformPrior(*this);
}
double UniformPrior::loglike(double x) const {
	return ll; // NB: does not check that x is between a and b... todo
}

BetaPrior::BetaPrior(double alpha, double beta) : alpha(alpha), beta(beta) {
	// TODO: check validity of params
	logBeta = gsl_sf_lnbeta(alpha, beta);
}
BetaPrior* BetaPrior::dup() const {
	return new BetaPrior(*this);
}
double BetaPrior::loglike(double x) const {
	return (alpha-1.0)*log(x) + (beta-1.0)*log(1.0-x) - logBeta;
}

NormalPrior::NormalPrior(double mu, double sigmasq) : mu(mu), sigmasq(sigmasq) {
	logsigmasq = log(sigmasq);
}
NormalPrior* NormalPrior::dup() const {
	return new NormalPrior(*this);
}
double NormalPrior::loglike(double x) const {
	return -0.5*(M_LN2 + M_LNPI + logsigmasq) - (x-mu)*(x-mu)/(2.0*sigmasq);
}

HalfNormalPrior::HalfNormalPrior(double sigmasq) : sigmasq(sigmasq) {
	logsigmasq = log(sigmasq);
}
HalfNormalPrior* HalfNormalPrior::dup() const {
	return new HalfNormalPrior(*this);
}
double HalfNormalPrior::loglike(double x) const {
	return -0.5*(M_LNPI - M_LN2 + logsigmasq) - x*x/(2.0*sigmasq);
}

StdNormalPrior::StdNormalPrior() { /* empty */ }
StdNormalPrior* StdNormalPrior::dup() const {
	return new StdNormalPrior();
}
double StdNormalPrior::loglike(double x) const {
	return -0.5*(M_LNPI + M_LN2) - x*x/2.0;
}

LogStdNormalPrior::LogStdNormalPrior() {
	/* empty */
}
LogStdNormalPrior* LogStdNormalPrior::dup() const {
	return new LogStdNormalPrior();
}
double LogStdNormalPrior::loglike(double x) const {
	double logx = log(x);
	return -0.5*(M_LNPI + M_LN2) - logx - logx*logx/2.0; // like is 1/(x * sqrt(2pi)) * exp(-log(x)^2/2)
}

GammaPrior::GammaPrior(double rate, double shape) : rate(rate), shape(shape) {
	// TODO: check validity of params
	lograte = log(rate);
	loggammashape = gsl_sf_lngamma(shape);
}
GammaPrior* GammaPrior::dup() const {
	return new GammaPrior(*this);
}
double GammaPrior::loglike(double x) const {
	return shape*lograte - loggammashape + (shape-1)*log(x) - rate*x;
}


InvGammaPrior::InvGammaPrior(double rate, double shape) : rate(rate), shape(shape) {
	// TODO: check validity of params
	lograte = log(rate);
	loggammashape = gsl_sf_lngamma(shape);
}
InvGammaPrior* InvGammaPrior::dup() const {
	return new InvGammaPrior(*this);
}
double InvGammaPrior::loglike(double x) const {
	return shape*lograte - loggammashape + (-shape-1)*log(x) - rate/x;
}

/* methods for transformations */

Transformation::Transformation() {
	/* empty */
}
Transformation::~Transformation() {
	/* empty */
}
double Transformation::evalFun(double x) const {
	return x;
}
double Transformation::evalJac(double x) const {
	return 1.0;
}
double Transformation::evalLogAbsJac(double x) const {
	return 0.0;
}


LinearTransformation::LinearTransformation(double loc, double scale) :
	loc(loc), scale(scale) {
	/* empty */
}
double LinearTransformation::evalFun(double x) const {
	return (x - loc) / scale;
}
double LinearTransformation::evalJac(double x) const {
	return 1.0/scale;
}
double LinearTransformation::evalLogAbsJac(double x) const {
	return -log(fabs(scale));
}


LogitTransformation::LogitTransformation(double loc, double scale) :
		loc(loc), scale(scale) {
	/* empty */
}
double LogitTransformation::evalFun(double x) const {
	return (log(x/(1.0-x)) - loc) / scale;
}
double LogitTransformation::evalJac(double x) const {
	return 1.0/(scale*x*(1.0-x));
}
double LogitTransformation::evalLogAbsJac(double x) const {
	return -log(fabs(scale*x*(1.0-x)));
}


LogTransformation::LogTransformation(double loc, double scale) :
		loc(loc), scale(scale) {
	/* empty */
}
double LogTransformation::evalFun(double x) const {
	return (log(x) - loc) / scale;
}
double LogTransformation::evalJac(double x) const {
	return 1/(x*scale);
}
double LogTransformation::evalLogAbsJac(double x) const {
	return -log(fabs(x*scale));
}
