#ifndef DISTRIBUTIONS_HPP
#define DISTRIBUTIONS_HPP

#include <utility> // pair

// likelihoods

double ll_binomial(int , double , int ); // n, p, k
double ll_poisson(double , int ); // lambda, k
double ll_neg_binomial(double , double , int ); // lambda=np/(1-p), n, k
double ll_polya(double , double , int ); // lambda=np/(1-p), d=1/n, k

double norm_lpdf(double x, double mu, double sigma);
double norm_lcdf(double x, double mu, double sigma);
double norm_lccdf(double x, double mu, double sigma);

std::pair<double, bool> censored_norm_lpdf(double x, double xhat, double sigma, int c);
std::pair<double, bool> censored_lognorm_lpdf(double x, double xhat, double sigma, int c);

std::pair<double, bool> binomial_lpmf(int k, int n, double p);
std::pair<double, bool> beta_binomial_lpmf(int k, int n, double p, double r);
std::pair<double, bool> poisson_lpmf(int k, double mu);
std::pair<double, bool> neg_binomial_lpmf(int k, double mu, double r);


// Priors

class Prior {
public:
	Prior();
	virtual ~Prior();
	virtual Prior* dup() const; // make a copy (borrowed from the D language)
	virtual double loglike(double ) const; // returns 0.0
protected:
};

class UniformPrior : public Prior {
public:
	UniformPrior(double , double );
	UniformPrior* dup() const override; // make a copy
	double loglike(double ) const;
protected:
	double a, b;
	double ll;
};

class BetaPrior : public Prior {
public:
	BetaPrior(double , double );
	BetaPrior* dup() const override; // make a copy
	double loglike(double ) const;
protected:
	double alpha, beta, logBeta;
};

class NormalPrior : public Prior {
public:
	NormalPrior(double , double ); // mu, sigmasq
	NormalPrior* dup() const override; // make a copy
	double loglike(double ) const;
protected:
	double mu, sigmasq, logsigmasq;
};

class StdNormalPrior : public Prior {
public:
	StdNormalPrior();
	StdNormalPrior* dup() const override;
	double loglike(double ) const;
protected:
};

class HalfNormalPrior : public Prior {
public:
	HalfNormalPrior(double ); // sigmasq
	HalfNormalPrior* dup() const override; // make a copy
	double loglike(double ) const;
protected:
	double sigmasq; // variance of the 'full' normal distribution
	double logsigmasq;
};

class LogStdNormalPrior : public Prior {
	LogStdNormalPrior();
	LogStdNormalPrior* dup() const override; // make a copy
	double loglike(double ) const;
protected:
};

class GammaPrior : public Prior {
public:
	GammaPrior(double , double ); // rate, shape
	GammaPrior* dup() const override;
	double loglike(double ) const;
protected:
	double rate, shape; // expectation is shape/rate
	double lograte, loggammashape; // pre-compute once!!
};

class InvGammaPrior : public Prior {
public:
	InvGammaPrior(double , double ); // rate, shape
	InvGammaPrior* dup() const override; // make a copy
	double loglike(double ) const;
protected:
	double rate, shape;
	double lograte, loggammashape; // pre-compute once!!
};

// Transformations Used for hypo-priors

class Transformation { // the base class acts as the identity
public:
	Transformation();
	virtual ~Transformation();
	virtual double evalFun(double ) const;
	virtual double evalJac(double ) const;
	virtual double evalLogAbsJac(double ) const;
protected:
};

class LinearTransformation : public Transformation {
public:
	LinearTransformation(double , double );
	double evalFun(double ) const override;
	double evalJac(double ) const override;
	double evalLogAbsJac(double ) const override;
protected:
	double loc, scale;
};

class LogitTransformation : public Transformation {
public:
	LogitTransformation(double , double );
	double evalFun(double ) const override;
	double evalJac(double ) const override;
	double evalLogAbsJac(double ) const override;
protected:
	double loc, scale;
};

class LogTransformation : public Transformation {
public:
	LogTransformation(double , double );
	double evalFun(double ) const override;
	double evalJac(double ) const override;
	double evalLogAbsJac(double ) const override;
protected:
	double loc, scale;
};

#endif
