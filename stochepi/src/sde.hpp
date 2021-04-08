#ifndef SDE_HPP_
#define SDE_HPP_

#include <boost/numeric/odeint.hpp>

#include "aux.hpp"

/** EM implementation taking std::vector.
 * Works with non-diagonal noise.
 */
class EulerMaruyamaStepper {
public:
  EulerMaruyamaStepper(Rng & rng) : rng(rng) { /* empty */ }
  typedef std::vector<double> state_type;
  typedef std::vector<double> deriv_type;
  typedef double value_type;
  typedef double time_type;
  typedef double order_type;
  typedef boost::numeric::odeint::stepper_tag stepper_category;
  static order_type order() {
    return 0.5;
  }
  template<class System>
  void do_step(System system, state_type & x, time_type t, time_type dt) const {
    deriv_type F(x.size()), G(x.size());
    system(x, F, t);
    // first query dimension of the Wiener process
    int wdim = system.get().wiener_dim();
    // sample a step of the Wiener process
    std::vector<double> Zi(wdim);
    std::generate(Zi.begin(), Zi.end(), [&](){return rng.normal(0, 1);});
    // now query for the dimension of the SDEs to allow for OD
    int cdim = system.get().continuous_dim();
    // sample a step of the Wiener process
    std::vector<double> Zo(cdim);
    std::generate(Zo.begin(), Zo.end(), [&](){return rng.normal(0, 1);});
    // use the volatility_vec_prod method from VectorField
    system.get().volatility_vec_prod(x, G, t, Zi, Zo);
    for ( size_t i = 0; i < x.size(); ++i ) {
      x[i] += dt * F[i] + sqrt(dt) * G[i]; // G = sigma * Z
    }
  }
private:
  /** ref to an RNG. Used to create normal deviates
   * @todo: construct actual RNG instead?
   */
  Rng & rng;
};

/** EM implementation taking std::vector.
 * Only works with diagonal noise.
 */
class EulerMaruyamaStepperDiag {
public:
  EulerMaruyamaStepperDiag(Rng & rng) : rng(rng) { /* empty */ }
  typedef std::vector<double> state_type;
  typedef std::vector<double> deriv_type;
  typedef double value_type;
  typedef double time_type;
  typedef double order_type;
  typedef boost::numeric::odeint::stepper_tag stepper_category;
  static order_type order() {
    return 0.5;
  }
  template<class System>
  void do_step(System system, state_type & x, time_type t, time_type dt) const {
    deriv_type F(x.size()), G(x.size());
    system(x, F, t);
    system.get().diffusion(x, G, t);
    for ( size_t i = 0; i < x.size(); ++i ) {
      double Z = rng.normal(0, 1);
      x[i] += dt * F[i] + sqrt(dt) * G[i] * Z;
    }
  }
private:
  /** ref to an RNG. Used to create normal deviates
   * @todo: construct actual RNG instead?
   */
  Rng & rng;
};


/** Runge-Kutta discretization of the Milstein scheme.
 * Only works for diagonal noise.
 *
 * @todo: check restrictions: we need that sigma_ii is
 * only dependent of x_i.
 */
 class RungeKuttaMilsteinStepper {
 public:
   RungeKuttaMilsteinStepper(Rng & rng) : rng(rng) { /* empty */ }
   typedef std::vector<double> state_type;
   typedef std::vector<double> deriv_type;
   typedef double value_type;
   typedef double time_type;
   typedef double order_type;
   typedef boost::numeric::odeint::stepper_tag stepper_category;
   static order_type order() {
     return 1.0;
   }
   /** The Millstein method takes steps as follows:
    * \f[
    *  Y_{n+1} = Y_n + f_n h + g_n \Delta W_n + \frac12 g_n g_n' [\Delta W_n^2 - h]
    * \f]
    * with \f$ f_n = f(Y_n) \f$, \f$ g_n = g(Y_n) \f$ and
    * \f$ g_n' = \frac{\partial g}{\partial x}(Y_n) \f$, and
    * \f$ \Delta W_n \sim \sqrt{h} \mathcal{N}(0,1) \f$.
    * The Kunge-Kutta approximation removes the need for a derivative,
    * and uses an intermediate step \f$ \bar{Y}_n = Y_n + f_n h + g_n \sqrt{h} \f$.
    * The next value is then calculated as
    * \f[
    *  Y_{n+1} = Y_n + f_n h + g_n \Delta W_n +
    *    \frac12 g_n \frac{1}{\sqrt{h}}(\bar{g}_n-g_n)[\Delta W_n^2 - h]
    * \f]
    * where \f$ \bar{g}_n = g(\bar{Y}_n) \f$.
    */
   template<class System>
   void do_step(System system, state_type & x, time_type t, time_type dt) const {
     deriv_type F(x.size()), G(x.size()), Gtilde(x.size());
     state_type xtilde(x.size());
     time_type sqdt = sqrt(dt);
     system(x, F, t);
     system.get().diffusion(x, G, t);
     // define intermediate step
     for ( size_t i = 0; i < x.size(); ++i ) {
       xtilde[i] = x[i] + dt*F[i] + sqdt*G[i];
     }
     // define Gtilde
     system.get().diffusion(xtilde, Gtilde, t);
     // update x
     for ( size_t i = 0; i < x.size(); ++i ) {
       double dW = rng.normal(0, 1) * sqdt;
       x[i] += dt*F[i] + G[i]*dW + 0.5*(Gtilde[i] - G[i]) * (dW*dW - dt) / sqdt;
     }
   }
 private:
   Rng & rng; // ref to an RNG
 };


/** @todo Implementation of an adaptive Stochastic Runge Kutta (SRK) method
 * using Rejection Sampling with Memory (RSwM)
 */





#endif
