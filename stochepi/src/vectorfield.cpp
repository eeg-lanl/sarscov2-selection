#include "vectorfield.hpp"

#include "variable.hpp"
#include "macros.hpp"



VectorField::VectorField(const State & s0, const std::vector<Transition*> & transitions,
    const Parameters & p) : p(p), s0(s0), dimension(0),
    continuous_dimension(0), wiener_dimension(0) {
  replace_transitions(transitions);
  // inspect transitions and group them into components of the vectorfield
  make_components();
  // assign indices to transitions and components in order to correctly integrate
  assign_indices();
}

VectorField::VectorField(const VectorField & vf) : p(vf.p), s0(vf.s0),
    dimension(0), continuous_dimension(0), wiener_dimension(0) { // copy constructor
  replace_transitions(vf.transitions);
  make_components(); // we have to re-make the components from the new transitions
  assign_indices(); // and then assign new indices
}

VectorField & VectorField::operator=(const VectorField & vf) { // copy assignment constructor
  if ( this != &vf ) {
    p = vf.p;
    s0 = vf.s0;
    replace_transitions(vf.transitions);
    make_components();
    assign_indices();
  } // else do nothing
  return *this;
}

void VectorField::make_components() {
  State zero = s0.zero_like();
  components.clear();
  components.resize(s0.flat_size());
  for ( auto trans : transitions ) {
    State epsilon = trans->apply(zero);
    for ( size_t i = 0; i < epsilon.flat_size(); ++i ) {
      if ( epsilon[i].disc_value != 0 ) {
        components[i].terms[trans] = epsilon[i].disc_value;
      }
    } // loop over State
  } // loop over transitions
}

void VectorField::assign_indices() {
  if ( s0.flat_size() != components.size() ) {
    throw std::logic_error("State and vectorfield have unequal dimension" + RIGHT_HERE);
  }
  for ( auto trans : transitions ) {
    // a priori, exclude transition from vector field
    trans->assignIndex(Indexable::invalid_index);
    // a priori, exclude transitions from the Wiener indices
    wiener_indices[trans] = Indexable::invalid_index;
  }
  // determine what vector field components to include in the vector field
  int idx = 0; // the index in the vector field
  int wiener_idx = 0;
  int continuous_idx = 0; // count the number of active components
  for ( size_t i = 0; i < s0.flat_size(); ++i ) {
    if ( !s0[i].is_discrete ) {
      components[i].assignIndex(idx++);
      /* terms that occur in one of the active components of the vector field
       * have to get a valid index in wiener_indices and contribute
       * to the wiener_dimension
       */
      for ( auto & term : components[i].terms ) {
        long & widx = wiener_indices.at(term.first); // NB: get a reference
        if ( widx == Indexable::invalid_index ) {
          widx = wiener_idx++; // This works because widx is a reference!
        }
        /* terms can occur in multiple components of the system of ODEs,
         * but we want to give it a wiener_idx only once!
         */
      } // for term in terms
      // count active components
      continuous_idx++;
    } else { // don't include component in the vector field
      components[i].assignIndex(Indexable::invalid_index);
      /* however, we need to include the transitions corresponding
       * to the component in the vector field to keep track of
       * stochasitic events
       */
      for ( auto & term : components[i].terms ) {
        if ( term.first->getIndex() == Indexable::invalid_index ) {
          term.first->assignIndex(idx++);
        }
        /* terms can occur in multiple components of the system of ODEs,
         * but we want to give it an index only once.
         */
      } // for loop over terms of component i
    } // else: discrete scalar
  } // for loop over the State-s variables
  // finally, set the dimension to the new value
  dimension = idx;
  continuous_dimension = continuous_idx;
  wiener_dimension = wiener_idx;
}

bool VectorField::update(const State & s0, const Parameters & p) {
  if ( this->s0.flat_size() != s0.flat_size() ) {
    throw std::invalid_argument("State is incompatible with vectorfield" + RIGHT_HERE);
  } // else... proceed without caution
  // re-define parameters
  this->p = p;
  // check that this s0 and s0 have different variable types
  bool update_indices = !compareSignatures(this->s0, s0);
  // re-define the initial state s0
  this->s0 = s0;
  // if necessary, re-assign indices
  if ( update_indices ) {
    assign_indices(); // assign indices also assignes the dimension
  }
  return update_indices;
}

int VectorField::dim() const {
  return dimension;
}

int VectorField::continuous_dim() const {
  return continuous_dimension;
}

int VectorField::discrete_dim() const {
  return dimension - continuous_dimension;
}

int VectorField::wiener_dim() const {
  return wiener_dimension;
}


void VectorField::operator()(const RealVec & y, RealVec & dydt,
    const double t) const {
  if ( y.size() != dydt.size() ) {
    throw std::invalid_argument(
      "state and total derivative do not have the same size" + RIGHT_HERE);
  }
  State s = decodeState(t, y);
  // 'first' part of vector: the actual derivatives of the variables
  for ( auto & comp : components ) {
    auto idx = comp.getIndex();
    if ( idx != Indexable::invalid_index ) {
      dydt[idx] = comp.rate(s, p);
    }
  }
  // 'second' part: derivatives of the transitions
  for ( auto trans : transitions ) {
    auto idx = trans->getIndex();
    if ( idx != Indexable::invalid_index ) {
      dydt[idx] = trans->rate(s, p); // make negative
    }
  }
  /* TODO when parameters are not time-dependent and all variables are
   * discrete, than we don't have to integrate anything!
   */
}

void VectorField::diffusion(const RealVec & y, RealVec & sig,
    const double t) const {
  if ( y.size() != sig.size() ) {
    throw std::invalid_argument(
      "state and total derivative do not have the same size" + RIGHT_HERE);
  }
  State s = decodeState(t, y);
  /* Diagonal noise on the components.
   * FIXME: use the off-diagonal elements of the diffusion matrix
   */
  // sig could not be initialized: set to zero!
  std::fill(sig.begin(), sig.end(), 0.0);
  for ( auto & comp : components ) {
    auto idx = comp.getIndex();
    if ( idx != Indexable::invalid_index ) {
      // diffusion term originating from discrete system approximation
      double diff = comp.diffusion(s, p);
      if ( diff < 0.0 ) {
        std::cerr << "# WARNING! diffusion term is negative: " << diff << std::endl;
        diff = 0.0;
      }
      // additional volatility term proportional to system size
      double volty = p[IDX_SIGMA_OD] * y[idx];
      // Assumption: diffusion and valatility are independent.
      sig[idx] = sqrt(diff + volty*volty);
    }
  }
}

void VectorField::volatility_vec_prod(const RealVec & y, RealVec & sig,
    const double t, const RealVec & Z) const {
  if ( y.size() != sig.size() ) {
    throw std::invalid_argument(
      "state and total derivative do not have the same size" + RIGHT_HERE);
  }
  if ( Z.size() != size_t(wiener_dim()) ) {
    throw std::invalid_argument(
      "Z must have length equal to the Wiener dimension" + RIGHT_HERE);
  }
  State s = decodeState(t, y);
  // sig could not be initialized: set to zero!
  std::fill(sig.begin(), sig.end(), 0.0);
  for ( auto & comp : components ) {
    auto idx = comp.getIndex();
    double elt = 0.0;
    if ( idx != Indexable::invalid_index ) {
      for ( auto & [trans, weight] : comp.terms ) {
        long widx = wiener_indices.at(trans);
        if ( widx == Indexable::invalid_index ) {
          throw std::logic_error("transition in an active vector field component has invalid Wiener index" + RIGHT_HERE);
        }
        double diff = trans->rate(s, p);
        if ( diff < 0 ) {
          std::cerr << "# WARNING! diffusion term is negative: " << diff << std::endl;
          diff = 0.0;
        }
        elt += weight * sqrt(diff) * Z[widx];
      } // for loop over terms of the vector field component
      sig[idx] = elt;
    }
  } // loop over components
}


void VectorField::volatility_vec_prod(const RealVec & y, RealVec & sig,
    const double t, const RealVec & Zi, const RealVec & Zo) const {
  if ( y.size() != sig.size() ) {
    throw std::invalid_argument(
      "state and total derivative do not have the same size" + RIGHT_HERE);
  }
  if ( Zi.size() != size_t(wiener_dim()) ) {
    throw std::invalid_argument(
      "Zi must have length equal to the Wiener dimension" + RIGHT_HERE);
  }
  if ( Zo.size() != size_t(continuous_dim()) ) {
    throw std::invalid_argument(
      "Zo must have length equal to the dimension" + RIGHT_HERE);
  }
  State s = decodeState(t, y);
  int cidx = 0;
  // sig could not be initialized: set to zero!
  std::fill(sig.begin(), sig.end(), 0.0);
  for ( size_t i = 0; i < components.size(); ++i ) {
    auto & comp = components[i];
    auto idx = comp.getIndex();
    double elt = 0.0;
    if ( idx != Indexable::invalid_index ) {
      for ( auto & [trans, weight] : comp.terms ) {
        long widx = wiener_indices.at(trans);
        if ( widx == Indexable::invalid_index ) {
          throw std::logic_error("transition in an active vector field component has invalid Wiener index" + RIGHT_HERE);
        }
        double diff = trans->rate(s, p);
        if ( diff < 0 ) {
          std::cerr << "# WARNING! diffusion term is negative: " << diff << std::endl;
          diff = 0.0;
        }
        elt += weight * sqrt(diff) * Zi[widx];
      } // for loop over terms of the vector field component
      // overdispersed noise
      double od_noise = 0.0;
      if ( !p.select_od.empty() ) { // use select_od to scale the noise
        od_noise = p[IDX_SIGMA_OD] * p.select_od[i] * y[idx] * Zo[cidx++];
      } else { // no selection vector included
        od_noise = p[IDX_SIGMA_OD] * y[idx] * Zo[cidx++];
      }
      sig[idx] = elt + od_noise;
    }
  } // loop over components
}



State VectorField::decodeState(double t, const RealVec & y) const {
  State s = s0; // make a copy of s0
  s.t() = t;
  for ( size_t i = 0; i < s0.flat_size(); ++i ) {
    // replace variables in the system of ODE with values from y
    auto idx = components[i].getIndex();
    if ( idx != Indexable::invalid_index ) {
      s[i] = y[idx];
    }
  }
  return s;
}


std::map<Transition*, double>
    VectorField::decodeLoads(double t, const RealVec & y) const {
  std::map<Transition*, double> loads;
  for ( auto trans : transitions ) {
    auto idx = trans->getIndex();
    if ( idx != Indexable::invalid_index ) {
      loads[trans] = y[idx];
    }
  }
  return loads;
}

double VectorField::totalTransitionRate() const {
  double lambda = 0.0;
  for ( auto trans : transitions ) {
    if ( trans->getIndex() != Indexable::invalid_index ) {
      lambda += trans->rate(s0, p);
    }
  }
  return lambda;
}



RealVec VectorField::encode(const State & s) const {
  RealVec y(dim(), 0.0); // TODO: encode remaining loads
  for ( size_t i = 0; i < components.size(); ++i ) {
    auto idx = components[i].getIndex();
    if ( idx != Indexable::invalid_index ) {
      y[idx] = s[i].value();
    }
  }
  return y;
}

// get drift and diffusion from the vector field

double VectorField::driftElt(const State & s, const Parameters & par, int i) const {
  if ( i < 0 || i >= int(components.size()) ) {
    throw std::range_error("index out of range" + RIGHT_HERE);
  } // else...
  return components[i].rate(s, par);
}

double VectorField::diffusionElt(const State & s, const Parameters & par, int i, int j) const {
  WARN_UNTESTED_FUN
  int n = components.size();
  if ( i < 0 || i >= n ) {
    throw std::range_error("index i out of range" + RIGHT_HERE);
  } // else...
  if ( j < 0 || j >= n ) {
    throw std::range_error("index j out of range" + RIGHT_HERE);
  } // else...
  if ( i == j ) {
    return components[i].diffusion(s, par);
  } // else, compute off-diagonal term
  double Dij = 0.0;
  for ( Transition* transition : transitions ) {
    auto it = components[i].terms.find(transition);
    if ( it != components[i].terms.end() ) {
      auto jt = components[j].terms.find(transition);
      if ( jt != components[j].terms.end() ) {
        Dij += transition->rate(s, par) * it->second * jt->second; // eta(s) * eps_i * eps_j
      }
    }
  }
  return Dij; // matrix element is zero
}





double VectorFieldComponent::rate(const State & s, const Parameters & par) const {
  double x = 0.0;
  for ( auto term : terms ) {
    // the second element of term is +/- 1 (in most cases), all rates are non-negative
    x += term.first->rate(s, par) * term.second;
  }
  return x;
}

double VectorFieldComponent::diffusion(const State & s, const Parameters & par) const {
  double x = 0.0;
  for ( auto term : terms ) {
    // sum absolute rates with squares of the coefficients of the increments
    x += term.first->rate(s, par) * term.second * term.second;
  }
  return x;
}
