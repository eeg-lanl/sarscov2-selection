#include "parprior.hpp"

#include <sstream>
#include <cmath> // fmod, sqrt...
#include <limits> // for infinity

#include "macros.hpp"

/* methods for ParPrior */

ParPrior::ParPrior() :
	name(""), elements(1), // use default value defined for ParElt
	lboundbool(false), lbound(0.0),
	uboundbool(false), ubound(0.0),
	homotopy_support(HomotopyClass::CONTRACTIBLE),
	prior(nullptr), rw_rule([](double t){return true;}),
	stepcounter(0), acceptcounter(0), ac(0), uc(1),
	selectionbool(false), selection(0),
	random_effects(false), loc(nullptr), scale(nullptr) {
	/* empty */
}

ParPrior::ParPrior(double value) :
	name(""), elements(1, ParElt(value)),
	lboundbool(false), lbound(0.0),
	uboundbool(false), ubound(0.0),
	homotopy_support(HomotopyClass::CONTRACTIBLE),
	prior(nullptr), rw_rule([](double t){return true;}),
	stepcounter(0), acceptcounter(0), ac(0), uc(1),
	selectionbool(false), selection(0),
	random_effects(false), loc(nullptr), scale(nullptr) {
	/* empty */
}

ParPrior::ParPrior(const std::vector<double> & values) :
	name(""), lboundbool(false), lbound(0.0),
	uboundbool(false), ubound(0.0),
	homotopy_support(HomotopyClass::CONTRACTIBLE),
	prior(nullptr), rw_rule([](double t){return true;}),
	stepcounter(0), acceptcounter(0), ac(0), uc(1),
	selectionbool(false), selection(0),
	random_effects(false), loc(nullptr), scale(nullptr) {
	// copy values to elements
	elements.resize(values.size());
	for ( size_t r = 0; r < values.size(); ++r ) {
		elements[r] = ParElt(values[r]);
	}
}

ParPrior::ParPrior(int n, double value) :
	name(""), elements(n, ParElt(value)),
	lboundbool(false), lbound(0.0),
	uboundbool(false), ubound(0.0),
	homotopy_support(HomotopyClass::CONTRACTIBLE),
	prior(nullptr), rw_rule([](double t){return true;}),
	stepcounter(0), acceptcounter(0), ac(0), uc(1),
	selectionbool(false), selection(0),
	random_effects(false), loc(nullptr), scale(nullptr) {
	/* empty */
}


// TODO: name!!
ParPrior::ParPrior(InitType ip, double value, double pstd) :
	name(""), elements(1, ParElt(value, pstd)),
	lboundbool(false), lbound(0.0),
	uboundbool(false), ubound(0.0),
	homotopy_support(HomotopyClass::CONTRACTIBLE),
	prior(nullptr), rw_rule([](double t){return true;}),
	stepcounter(0), acceptcounter(0), ac(0), uc(1),
	selectionbool(false), selection(0),
	random_effects(false), loc(nullptr), scale(nullptr) { // some common constructions
	switch ( ip ) {
		case UNRESTRICTED: {
			break;
		}
		case POSITIVE: {
			lboundbool = true;
			lbound = 0.0;
			break;
		}
		case NEGATIVE: {
			uboundbool = true;
			ubound = 0.0;
			break;
		}
		default: {
			throw std::logic_error("invalid InitType" + RIGHT_HERE);
			break; // redundant
		}
	}
}


ParPrior::ParPrior(double value, double pstd, Prior* prior, std::string name) :
	name(name), elements(1, ParElt(value, pstd)),
	lboundbool(false), lbound(0.0),
	uboundbool(false), ubound(0.0),
	homotopy_support(HomotopyClass::CONTRACTIBLE),
	prior(prior), rw_rule([](double t){return true;}),
	stepcounter(0), acceptcounter(0), ac(0), uc(1),
	selectionbool(false), selection(0),
	random_effects(false), loc(nullptr), scale(nullptr) {
	/* empty */
}

ParPrior::ParPrior(const ParPrior & pp) {
	/* NB: don't call clear(), since members are not initialized.
	 * in particular: pointers are not NULL
	 */
	copy(pp);
}

ParPrior::~ParPrior() {
	clear();
}

ParPrior & ParPrior::operator=(const ParPrior & pp) {
	if ( this != &pp ) {
		clear(); // deletes pointers to prior, loc ans scale
		copy(pp); // copies data from pp
	}
	return *this;
}

void ParPrior::copy(const ParPrior & pp) {
	name = pp.name;
	elements = pp.elements;
	lboundbool = pp.lboundbool;
	lbound = pp.lbound;
	uboundbool = pp.uboundbool;
	ubound = pp.ubound;
	stepcounter = pp.stepcounter;
	acceptcounter = pp.acceptcounter;
	ac = pp.ac;
	uc = pp.uc;
	if ( pp.prior != nullptr ) {
		prior = pp.prior->dup();
	}	else {
		prior = nullptr;
	}
	homotopy_support = pp.homotopy_support;
	rw_rule = pp.rw_rule;
	selectionbool = pp.selectionbool;
	selection = pp.selection;
	// copy random effects data
	random_effects = pp.random_effects;
	if ( pp.random_effects ) {
		// pp.random_effects == true guarentees that pp.loc isn't NULL or garbage
		loc = pp.loc->dup();
		scale = pp.scale->dup();
	}	else {
		loc = nullptr;
		scale = nullptr;
	}
}

void ParPrior::clear() {
	delete prior;
	delete loc;
	delete scale;
}


ParPrior & ParPrior::operator=(double new_value) {
	if ( !inBounds(new_value) ) {
		throw std::invalid_argument("new value is not within range" + RIGHT_HERE);
	}
	for ( auto & elt : elements ) {
		elt.old_value = elt.value;
		elt.value = new_value;
	}
	return *this;
}

ParPrior & ParPrior::operator=(const std::vector<double> & new_values) {
	elements.resize(new_values.size());
	for ( size_t r = 0; r < new_values.size(); ++r ) {
		if ( !inBounds(new_values[r]) ) {
			throw std::invalid_argument("new value is not within range" + RIGHT_HERE);
		}
		elements[r] = ParElt(new_values[r]);
	}
	// turn off previous selection, as it may no longer be valid
	selectionbool = false;
	selection = 0;
	return *this;
}

ParPrior & ParPrior::operator*=(double x) {
	for ( auto & elt : elements ) {
		double new_value = elt.value * x;
		if ( !inBounds(new_value) ) {
			throw std::invalid_argument("new value is not within range" + RIGHT_HERE);
		}
		elt.old_value = elt.value;
		elt.value = new_value;
	}
	return *this;
}

ParPrior & ParPrior::operator/=(double x) {
	if ( x == 0 ) {
		throw std::invalid_argument("denominator is zero" + RIGHT_HERE);
	}
	for ( auto & elt : elements ) {
		double new_value = elt.value / x;
		if ( !inBounds(new_value) ) {
			throw std::invalid_argument("new value is not within range" + RIGHT_HERE);
		}
		elt.old_value = elt.value;
		elt.value = new_value;
	}
	return *this;
}

ParPrior & ParPrior::operator+=(double x) {
	for ( auto & elt : elements ) {
		double new_value = elt.value + x;
		if ( !inBounds(new_value) ) {
			throw std::invalid_argument("new value is not within range" + RIGHT_HERE);
		}
		elt.old_value = elt.value;
		elt.value = new_value;
	}
	return *this;
}

ParPrior & ParPrior::operator-=(double x) {
	for ( auto & elt : elements ) {
		double new_value = elt.value - x;
		if ( !inBounds(new_value) ) {
			throw std::invalid_argument("new value is not within range" + RIGHT_HERE);
		}
		elt.old_value = elt.value;
		elt.value = new_value;
	}
	return *this;
}



/** Use the getValue method to retrieve the correct value
 * which can depend on isLocked()
 */
ParPrior::operator double() const {
	return getValue();
}

/** Manually get one of the values */
const double & ParPrior::operator[](size_t r) const {
	if ( r >= elements.size() ) {
		throw std::range_error("selection not in range of elements" + RIGHT_HERE);
	} // else...
	return elements[r].value;
}

double & ParPrior::operator[](size_t r) {
	if ( r >= elements.size() ) {
		throw std::range_error("selection not in range of elements" + RIGHT_HERE);
	} // else...
	return elements[r].value;
}



ParPrior* ParPrior::dup() const {
	return new ParPrior(*this);
}

void ParPrior::lock() {
	// keep locked tot current value
	for ( auto & elt : elements ) {
		elt.locked = true;
	}
	if ( random_effects ) {
		loc->lock(); scale->lock();
	}
}

void ParPrior::lockAllBut(int s) {
	// locks all elements to current value, unlocks r-th element
	for ( int r = 0; r < int(elements.size()); ++r ) {
		elements[r].locked = (r != s);
	}
}

void ParPrior::removeSingleLocks() {
	if ( !isLocked() ) {
		unlock();
	}
}


void ParPrior::unlock() {
	for ( auto & elt : elements ) {
		elt.locked = false;
	}
	if ( random_effects ) {
		loc->unlock(); scale->unlock();
	}
}

bool ParPrior::isLocked() const {
	// are all parameters locked?
	for ( auto & elt : elements ) {
		if ( !elt.locked ) {
			return false; // return false if ANY one of the parameters is not locked
		}
	}
	if ( random_effects && (!loc->isLocked() || !scale->isLocked()) ) {
		// NB: && operator short-circuits
		return false;
	}
	// at reaching this point, we know that all elements are locked
	return true;
}

void ParPrior::setPrior(Prior* new_prior) {
	delete prior;
	prior = new_prior;
}




void ParPrior::setBounds(double lbound, double ubound) {
	this->lbound = lbound;
	lboundbool = true;
	this->ubound = ubound;
	uboundbool = true;
	// check validity
	for ( auto & elt : elements ) {
		if ( !inBounds(elt.value) ) {
			std::stringstream ss;
			ss << "given bounds for '" << getName() << "' do not contain value " << elt.value << RIGHT_HERE;
			throw std::logic_error(ss.str());
		}
	}
}

void ParPrior::setNameBoundsPstdAndUnlock(std::string name, double lbound,
		double ubound, double pstd) {
	setName(name);
	setBounds(lbound, ubound);
	setPstd(pstd);
	unlock();
}

void ParPrior::setRWRule(RWRule rw_rule) {
	this->rw_rule = rw_rule;
}

void ParPrior::setLBound(double lbound) {
	this->lbound = lbound;
	lboundbool = true;
	// check validity
	for ( auto & elt : elements ) {
		if ( !inBounds(elt.value) ) {
			std::stringstream ss;
			ss << "given bounds for '" << getName() << "' do not contain value " << elt.value << RIGHT_HERE;
			throw std::logic_error(ss.str());
		}
	}
}

void ParPrior::setUBound(double ubound) {
	this->ubound = ubound;
	uboundbool = true;
	// check validity
	for ( auto & elt : elements ) {
		if ( !inBounds(elt.value) ) {
			std::stringstream ss;
			ss << "given bounds for '" << getName() << "' do not contain value " << elt.value << RIGHT_HERE;
			throw std::logic_error(ss.str());
		}
	}
}

void ParPrior::setRandomEffects(double loc_val, double loc_pstd, double scale_val, double scale_pstd) {
	if ( random_effects ) {
		// first delete any old loc and scale
		delete loc; delete scale;
	}
	random_effects = true;
	loc = new ParPrior(UNRESTRICTED, loc_val, loc_pstd);
	loc->setName(name + "__loc");
	scale = new ParPrior(POSITIVE, scale_val, scale_pstd);
	scale->setName(name + "__scale");
	setPrior(new StdNormalPrior()); // TODO: other options...
}

bool ParPrior::isRandomEffects() const {
	return random_effects;
}

const ParPrior* ParPrior::getLoc() const {
	return loc;
}
const ParPrior* ParPrior::getScale() const {
	return scale;
}

ParPrior* ParPrior::getLoc() {
	return loc;
}
ParPrior* ParPrior::getScale() {
	return scale;
}



bool ParPrior::isBounded() const {
  return lboundbool && uboundbool;
}

bool ParPrior::inBounds(double x) const {
	// check if the argument is between lbound and ubound
	return ((!uboundbool || x <= ubound) && (!lboundbool || x >= lbound));
}

double ParPrior::getLengthInterval() const {
	if ( isBounded() ) {
    return ubound - lbound;
  } else {
    return std::numeric_limits<double>::infinity();
  }
}

void ParPrior::setHomotopy(HomotopyClass hc) {
	switch ( hc ) {
		case HomotopyClass::CONTRACTIBLE: {
			homotopy_support = HomotopyClass::CONTRACTIBLE;
			break;
		}
		case HomotopyClass::CIRCULAR: {
			if ( lboundbool && uboundbool ) {
				homotopy_support = HomotopyClass::CIRCULAR;
			}
			else {
				throw std::logic_error("can't make an unbounded support CIRCULAR" + RIGHT_HERE);
			}
			break;
		}
		default: {
			throw std::logic_error("invalid HomotopyClass given" + RIGHT_HERE);
			break;
		}
	}
}


void ParPrior::select(int r) {
	if ( r < 0 || r >= int(elements.size()) ) {
		throw std::range_error("selection not in range of elements" + RIGHT_HERE);
	} // else...
	selectionbool = true;
	selection = r;
}

void ParPrior::deselect() {
	selectionbool = false;
	selection = 0;
}


double ParPrior::getValue() const {
	if ( size() > 1 ) {
		if ( selectionbool ) {
			return elements[selection].value;
		} else {
			std::cout << *this << std::endl; // TESTING!
			throw std::logic_error("call to getValue() of vector-valued ParPrior is ambiguous" + RIGHT_HERE);
		}
	} else { // scalar
		return elements[0].value;
	}
}

void ParPrior::setPstd(double pstd) {
	for ( auto & elt : elements ) {
		elt.pstd = pstd;
	}
}

double ParPrior::getAr() const {
	if ( stepcounter > 0 ) {
		return double(acceptcounter) / stepcounter;
	}	else {
		return 0.0;
	}
}


double ParPrior::loglike() const {
	if ( random_effects ) {
		LinearTransformation fun(loc->getValue(), scale->getValue());
		return loglike(&fun);
	} else {
		return loglike(nullptr);
	}
}


double ParPrior::loglike(Transformation* fun) const {
	double ll = 0.0;
	for ( auto & elt : elements ) {
		if ( prior != nullptr ) {
			if ( fun == nullptr ) {
				ll += prior->loglike(elt.value);
			}	else {
				ll += prior->loglike(fun->evalFun(elt.value)) + fun->evalLogAbsJac(elt.value);
			}
		}	else if ( fun != nullptr ) {
			// assume that prior is not informative
			ll += fun->evalLogAbsJac(elt.value);
		}
	}
	if ( random_effects ) {
		ll += loc->loglike() + scale->loglike();
	}
	return ll;
}

/** if lboundbool, mirror in lbound. if uboundbool, mirror in ubound
 * if both uboundbool and lboundbool, mirror in both bounds,
 * by using a combination of mirroring and modulus
 * the optional parameter rel_temp can be used for cooling-down schemes
 * e.g. in simulated annealing
 */
void ParPrior::mutate(Rng & rng, double rel_temp, double t) {
	// first check if the rw_rule allows mutation...
	if ( !rw_rule(t) ) return;
	// then run through individual elements
	for ( auto & elt : elements ) {
		if ( elt.locked ) continue;
		// else...
		elt.old_value = elt.value;
		elt.value = elt.value + rng.stdNormal() * elt.pstd * rel_temp; // normal proposal
		if ( lboundbool ) {
			if ( uboundbool ) { // the 'difficult' case
				double interval = (ubound-lbound); // length interval
				switch ( homotopy_support ) {
					case HomotopyClass::CONTRACTIBLE: { // apply double mirroring
						elt.value = fmod(elt.value-lbound, 2*interval);
						// remainder after dividion by denominator
						if ( elt.value < 0.0 ) elt.value += 2*interval;
						// fmod does not always return non-negative numbers
						if ( elt.value > interval ) elt.value = 2*interval - elt.value;
						// mirror
						elt.value += lbound;
						// translate back to the 'true' interval
						break;
					}
					case HomotopyClass::CIRCULAR: { // value modulo interval
						elt.value = fmod(elt.value-lbound, interval);
						if ( elt.value < 0.0 ) elt.value += interval;
						// fmod can return negative
						elt.value += lbound; // translate back
						break;
					}
					default: { // POINT handled above...
						throw std::logic_error("invalid HomotopyClass" + RIGHT_HERE);
						break;
					}
				}
			}	else { // mirror in lbound?
				if ( elt.value < lbound ) {
					elt.value = 2*lbound - elt.value;
				}
			}
		}	else { // lboundbool == false!!
			if ( uboundbool ) { // mirror in ubound?
				if ( elt.value > ubound ) {
					elt.value = 2*ubound - elt.value;
				}
			}
			// else: uboundbool == false && uboundbool == false, so value is valid!
		}
	} // for elt in elements
	if ( random_effects ) {
		loc->mutate(rng, rel_temp, t); scale->mutate(rng, rel_temp, t);
	}
	stepcounter++;
}

void ParPrior::accept() {
	if ( !isLocked() ) {
		ac++; // used by updatePvar
		acceptcounter++; // used for diagnostics
	}
	if ( random_effects ) {
		loc->accept(); scale->accept();
	}
}

void ParPrior::reject() {
	for ( auto & elt : elements ) {
		if ( !elt.locked ) {
			elt.value = elt.old_value;
		}
	}
	if ( random_effects ) {
		loc->reject(); scale->reject();
	}
}

void ParPrior::updatePstd() {
	for ( auto & elt : elements ) {
		if ( !elt.locked ) {
			if ( stepcounter % PVAR_UPDATE_INTERVAL == 0 ) {
				double ar = double(ac) / PVAR_UPDATE_INTERVAL;
				// lower, or increase the proposal variance
				if ( ar < OPTIMAL_ACCEPTANCE_RATE ) {
					elt.pstd *= exp(-0.5/sqrt(uc)); // jumps are too big!
				}	else {
					elt.pstd *= exp(0.5/sqrt(uc)); // jumps are too small!
				}
				// do some checks to prevent explosions
				if ( isBounded() ) {
					elt.pstd = std::min(elt.pstd, getLengthInterval());
				}
				// update/reset counters
				ac = 0; // reset acceptance counter
				uc++; // increase update counter
			}
		} // if ! locked... else do nothing
	}
	if ( random_effects ) {
		loc->updatePstd(); scale->updatePstd();
	}
}

std::string ParPrior::getName() const {
	return name;
}

void ParPrior::setName(std::string name) {
	this->name = name;
}

size_t ParPrior::size() const {
	return elements.size();
}

void ParPrior::print(std::ostream & os) const {
	os << "<param " << "name='" << name << "' "
	   << "re='" << std::boolalpha << random_effects << std::noboolalpha << "' "
	   << ">" << std::endl;
	for ( size_t r = 0; r < elements.size(); ++r ) {
		auto & elt = elements[r]; // alias
		os << "<elt "
		   << "idx='" << r << "' "
		   << "val='" << elt.value << "' "
		 	 << "pstd='" << elt.pstd << "' "
		 	 << "lock='" << std::boolalpha << elt.locked << std::noboolalpha << "' "
			 << "/>" << std::endl;
	}
	// don't print hyper parameters: taken care of by Parameters class
	os << "</param>";
}




// non menber-functions for ParPrior

std::istream & operator>>(std::istream & is, ParPrior & par) {
  // read value from stream
	double value;
	is >> value; // use standard stream for double
	par = value; // use ParPrior::operator=(double )
  return is;
}
