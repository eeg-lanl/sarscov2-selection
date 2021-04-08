#ifndef VECTORFIELD_HPP_
#define VECTORFIELD_HPP_

#include <map> // for map Transition* -> int (weight)
#include <vector>
#include <unordered_map> // for map Transition* -> size_t (Wiener process index)

#include "state.hpp"
#include "parameter.hpp"
#include "transition.hpp"
#include "aux.hpp"


class VectorFieldComponent : public Indexable {
public:
  double rate(const State & , const Parameters & ) const;
  // the diffusion term (diagonal) in the Fokker-Planck equation
  double diffusion(const State & , const Parameters & ) const;
  friend class VectorField; // VectorField needs access to terms
protected:
  std::map<Transition*, int> terms;
};

class VectorField : public TransitionList {
public:
  VectorField(const State & , const std::vector<Transition*> & ,
      const Parameters & );
  VectorField(const VectorField & ); // copy constructor
  VectorField & operator=(const VectorField & ); // copy assignment constructor
  /* after switch discrete <--> continuous, we have to re-assign
   * indices. If indices were re-assigned, "true" is returned, else "false".
   * the dimension can change after a call to update.
   * Update also re-defines parameters.
   */
  bool update(const State & , const Parameters & );
  int dim() const; // number of ODEs, including cumulative rates
  int continuous_dim() const; // dimension of the actual vector field: the SDEs or ODEs
  int discrete_dim() const; // number of discrete events
  int wiener_dim() const; // required dimension of the Wiener process
  /** define the operator() in order to use boost odeint.
   * This is also the drift term if the VectorField encodes a system
   * of SDEs.
   */
  void operator()(const std::vector<double> & y,
      std::vector<double> & dydt, const double t) const;
  /** Diffusion term in system of SDEs */
  void diffusion(const std::vector<double> & y,
      std::vector<double> & dydt, const double t) const;
  /** Used for non-diagonal diffusion in SDEs */
  void volatility_vec_prod(const RealVec & y, RealVec & sig,
      const double t, const RealVec & Z) const;
  /** Used for non-diagonal diffusion in SDEs,
   * now allowing for overdispersed noise
   */
  void volatility_vec_prod(const RealVec & y, RealVec & sig,
      const double t, const RealVec & Zi, const RealVec & Zo) const;
  // @todo: methods for the cumulative/instantaneous transition rate
  // use VectorField to retrieve the final state and transition loads after integration
  State decodeState(double t, const RealVec & ) const;
  std::map<Transition*, double> decodeLoads(double t, const RealVec & ) const;
  /** get the sum of the transition rates in order to anticipate a good
   * time step. Takes s0 as state and p as parameters.
   */
  double totalTransitionRate() const;
  std::vector<double> encode(const State & s) const; // for initial conditions
  // elements of the diffusion matrix and drift vector
  double driftElt(const State &, const Parameters &, int i) const;
  /** Get the (i,j)-th component of the diffusion matrix G.
   * The Fokker-Planck equation is given by
   * \f[
   *  \frac{\partial u}{\partial t} = -\frac{\partial}{\partial x}[F u] +
   *  \frac12 \frac{\partial^2}{\partial x^2}[Gu]
   * \f]
   */
  double diffusionElt(const State &, const Parameters &, int i, int j) const;
protected:
  // auxiliary functions
  void make_components();
  // assign_indices also sets the dimension (as returned by dim())
  void assign_indices();
  std::vector<VectorFieldComponent> components;
  std::unordered_map<Transition*, long> wiener_indices;
  Parameters p; // p is required for call to operator()
  State s0; // the initial state
  int dimension; // set by assign_indices
  int continuous_dimension; // set by assign_indices: number of active components
  int wiener_dimension; // set by assign_indices
};



#endif
