#ifndef PARPRIOR_HPP
#define PARPRIOR_HPP

#include <iostream>
#include <string>
#include <vector> // for vector-valued parameters
#include <functional>

#include "aux.hpp" // interfaces
#include "rng.hpp"
#include "distributions.hpp" // Prior


/** here we define objects that serve as parameters of the models.
 * In the current algorithm, the parameter class has a pointer
 * to a prior distribution, and handles the proposal distribution.
 * TODO: 1) make one class for discrete and continuous variables\
 * 2) allow for parameter vectors and multivariate proposals.
 */

constexpr int PVAR_UPDATE_INTERVAL = 50; // TODO: parameter
constexpr double OPTIMAL_ACCEPTANCE_RATE = 0.44; // 0.44 for Metropolis within Gibbs with normal proposal


enum class HomotopyClass {
	CONTRACTIBLE=0,
	CIRCULAR,
};

/** Rule based on time to lock parameters during certain parts
 * of the timeseries.
 * @todo: add Parameters as a second argument.
 */
typedef std::function<bool(double)> RWRule;


class ParPrior : public Printable {
public:
	enum InitType {
		UNRESTRICTED,
		POSITIVE,
		NEGATIVE
	};
	ParPrior(); // constant 0.0 and some defaults
	ParPrior(double ); // sets value to argument, other members to defaults
	ParPrior(const std::vector<double> & );
	ParPrior(int , double );
	ParPrior(InitType , double , double ); // some common constructions
	ParPrior(double , double , Prior* , std::string name="");
	// initial value, proposal variance, a prior, and an optional name
	ParPrior(const ParPrior & );
	~ParPrior(); // deletes prior, loc and scale
	// copy of prior, uses Prior::dup method
	ParPrior & operator=(const ParPrior & );
	// replace the value with something new
	ParPrior & operator=(double );
	ParPrior & operator=(const std::vector<double> & );
	// some arithmatic...
	ParPrior & operator*=(double );
	ParPrior & operator/=(double );
	ParPrior & operator+=(double );
	ParPrior & operator-=(double );
	// conversion to double
	operator double() const;
	// get element
	const double & operator[](size_t r) const;
	double & operator[](size_t r);
	// copy of prior, uses Prior::dup method
	virtual ParPrior* dup() const; // make a copy and return pointer
	void lock(); // keep locked to current value
	void lockAllBut(int r); // locks all elements, unlocks r-th element
	void removeSingleLocks(); // IF some elements are not locked, remove ALL locks
	void unlock(); // remove the lock
	/** convenient combination of the functions setName(), setBounds(),
	 * setPstd() and unlock(), which often occur together in the code.
	 */
	void setNameBoundsPstdAndUnlock(std::string name, double lbound,
			double ubound, double pstd);
	void setRWRule(RWRule rw_rule);
	bool isLocked() const; // false is ANY element is unlocked
	void setPrior(Prior* new_prior);
	void setBounds(double , double );
	void setLBound(double );
	void setUBound(double );
	void setRandomEffects(double loc_val, double loc_pstd,
			double scale_val, double scale_pstd);
	bool isRandomEffects() const;
	const ParPrior* getLoc() const;
	const ParPrior* getScale() const;
	ParPrior* getLoc();
	ParPrior* getScale();
	void setPstd(double );
	bool isBounded() const; // uboundbool && lboundbool
	bool inBounds(double ) const; // check if the argument is between lbound and ubound
	double getLengthInterval() const; // inf if unbounded...
	void setHomotopy(HomotopyClass );
	double getValue() const;
	void select(int r);
	void deselect();
	double getAr() const; // total acceptance ratio
	double loglike() const;
	double loglike(Transformation* fun) const;
	void mutate(Rng & rng, double rel_temp, double t);
	void accept(); // increase the acceptance counter
	void reject(); // replace value with old_value
	void updatePstd();
	void resetUpdateCounter(); // sets uc to 0: allows pvar to be re-optimized.
	std::string getName() const;
	void setName(std::string name);
	size_t size() const;
	void print(std::ostream & ) const override; // implement Printable interface
protected:
	// internal struct for basic parameter object
	struct ParElt {
		ParElt() : value(0.0), old_value(0.0), pstd(1.0), locked(true) { /* empty */ }
		ParElt(double x) : value(x), old_value(x), pstd(1.0), locked(true) { /* empty */ }
		ParElt(double x, double s) : value(x), old_value(x), pstd(s), locked(true) { /* empty */ }
		double value;
		double old_value;
		double pstd;
		bool locked;
	};
	// data
	std::string name;
	std::vector<ParElt> elements;
	bool lboundbool;
	double lbound;
	bool uboundbool;
	double ubound;
	// TODO: implement Stan-like transformations (log or logit) to achieve bounded parameters
	HomotopyClass homotopy_support; // contractible, circular, ...
	Prior* prior; // the prior distribution for the elements
	RWRule rw_rule;
	int stepcounter; // counts number of steps taken
	int acceptcounter; // overall
	int ac; // acceptance counter (short time scale)
	int uc; // update counter
	// index of the value to be returned
	bool selectionbool;
	int selection;
	/* the vector-valued parameter object can be used to model random effects
	 * in this case the loc and scale members are used to compute the likelihood
	 */
	bool random_effects; // random effects guards loc and scale: they have been allocated ONLY IF random_effects is true
	ParPrior *loc, *scale;
	// auxiliary methods
	void copy(const ParPrior & pp); // for copy constructor and copy assignment constructor
	void clear(); // deletes pointers
};

// istream operator for ParPrior
std::istream & operator>>(std::istream & is, ParPrior & par);

#endif
