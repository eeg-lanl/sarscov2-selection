#ifndef _PARAMETER_HPP
#define _PARAMETER_HPP

#include <vector>
#include <map>
#include <list>

#include "aux.hpp" // Printable
#include "parprior.hpp"
#include "parse_params.hpp" // ParSpecs

// HACK 1: the first parameter is always the overdispersion of the noise
constexpr size_t IDX_SIGMA_OD = 0;

class Parameters : public Printable {
public:
  Parameters() { /* empty */ }
  Parameters(int n);
  Parameters(const std::vector<size_t> & shape);
  // TODO initilizer list
  const ParPrior & operator[](int i) const; // retrieve a parameter by it's index (can use enum)
  ParPrior & operator[](int i); // retrieve a parameter by it's index (can use enum)
  const ParPrior & operator[](std::string name) const; // retrieve a parameter by it's name
  ParPrior & operator[](std::string name); // retrieve a parameter by it's name
  const ParPrior & operator()(size_t i) const; // retrieve i-th parameter, must be scalar
  ParPrior & operator()(size_t i); // retrieve i-th parameter, must be scalar
  const ParPrior & operator()(size_t i, size_t j) const; // retrieve (i,j)-th parameter according to shape
  ParPrior & operator()(size_t i, size_t j); // retrieve (i,j)-th parameter according to shape
  bool isScalar(size_t i) const; // is the i-th parameter a scalar according to shape?
  size_t flat_size() const; // returns paramvec.size()
  size_t size() const; // returns shape.size()
  void mutate(Rng & rng, double rel_temp, double t);
  void select(int r); // make sure that operator() returns the r-th element of vector-valued parameters
  void deselect(); // remove selection
  void lockAllBut(int r); // locks all parameter elements except the r-th if at least one element is NOT locked
  void removeSingleLocks(); // removes all locks for random effects parameters if at least one element is NOT locked
  void print(std::ostream & os) const override;
  double loglike() const; // prior likelihood
  // HACK 2: we have to disable some of the overdispersed noise terms
  std::vector<double> select_od;
protected:
  // aux methods for shape
  std::pair<size_t, size_t> unflatten_index(size_t i) const;
  size_t flatten_index(size_t i, size_t j) const;
  // protected members
  std::vector<ParPrior> paramvec;
  std::vector<size_t> shape;
};


template <class T>
void load_param_specs(Parameters & par, const std::map<T, std::string> & parNameMap,
      const std::list<ParSpecs> & parSpecsList) {
  for ( auto & [symb, name] : parNameMap ) {
    // find parameter in parSpecs
    auto it = std::find_if(parSpecsList.begin(), parSpecsList.end(),
        [name](const ParSpecs & psc) -> bool {return psc.name == name;});
    if ( it == parSpecsList.end() ) {
      throw std::runtime_error("parameter '" + name +
        "' is not defined in the parameter specs list" + RIGHT_HERE);
    } // else... parameter is defined
    ParSpecs psc = *it;
    // set value
    par[symb] = psc.value;
    // set name
    par[symb].setName(name);
    // if not locked, set bounds and sigma for the random walk
    if ( !psc.is_locked ) {
      par[symb].unlock();
      par[symb].setPstd(psc.sigma_rw);
      if ( psc.is_lbound ) {
        par[symb].setLBound(psc.lbound);
      }
      if ( psc.is_ubound ) {
        par[symb].setUBound(psc.ubound);
      }
    }
  }
}



#endif
