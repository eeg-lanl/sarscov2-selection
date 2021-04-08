#ifndef SARS2MUT_HPP_
#define SARS2MUT_HPP_

#include "transition.hpp"

namespace sarsmodel {
  enum VarNames {
    S=0,
    Ew,
    Em,
    Iw,
    Im,
    H, // hospitalizations (or severe infections)
    R, // removed / recovered (or dead)
    dIw, // auxilliary cumulative incidence for case incidence data
    dIm, // auxilliary cumulative incidence for case incidence data
    dD, // auxilliary cumulative incidence for death incidence data
    NUM_VARIABLES
  };

  enum ParSymbol {
    // index 0 for overdispersion parameter
    t0=1, // initial time
    N0, // initial population size
    lambda, // import rate of infected individuals
    q_mut, // fraction of mutant imports
    p_mut, // fraction of initial mutant fraction
    sigma, // fitness advantage of the mutant
    beta0, // infection rate 0
    beta1, // infection rate 1
    beta2, // infection rate 2
    beta3, // infection rate 3
    t1, // first breakpoint
    t2, // snd breakpoint
    t3, // trd breakpoint
    upsilon1, // time window to transition between conditions 0 and 1
    upsilon2, // time window to transition between conditions 1 and 2
    upsilon3, // time window to transition between conditions 2 and 3
    alpha, // transition from exposed to infectious
    gamma, // recovery rate
    nu, // hospitalization rate / severe infection
    omega, // death or recovery from severe infection
    delta, // probability of death after severe infection
    mu, // mutation rate
    theta, // testing rate
    theta_r, // overdispersion of reporting
    eta_r, // overdispersion of sequencing
    delta_r, // overdispersion of deaths
    epsilon, // determines initial condition
    xi, // fraction inititally immune/recovered
    NUM_PARAMETERS
  };

  extern const std::map<ParSymbol, std::string> parNames;

  State observables(const State & s, const Parameters & par);

  double infection_rate(const State & s, const Parameters & par);

  double testing_rate(const State & s, const Parameters & par);

  // HACK: make sure that we can replace the import freq with a spline
  extern std::function<double(double, const Parameters & par)> import_freq_mut;
  extern std::function<double(double, const Parameters & par)> import_rate;

  State gen_init_state(const Parameters &, Rng & rng);

  // transitions

  class TransSEw : public CloneableTransition<TransSEw> {
  public:
    TransSEw() : CloneableTransition("TransSEw") {}
    State apply(const State & ) const override;
    double rate(const State & , const Parameters & ) const override;
  };

  class TransSEm : public CloneableTransition<TransSEm> {
  public:
    TransSEm() : CloneableTransition("TransSEm") {}
    State apply(const State & ) const override;
    double rate(const State & , const Parameters & ) const override;
  };

  class TransEwIw : public CloneableTransition<TransEwIw> {
  public:
    TransEwIw() : CloneableTransition("TransEwIw") {}
    State apply(const State & ) const override;
    double rate(const State & , const Parameters & ) const override;
  };

  class TransEmIm : public CloneableTransition<TransEmIm> {
  public:
    TransEmIm() : CloneableTransition("TransEmIm") {}
    State apply(const State & ) const override;
    double rate(const State & , const Parameters & ) const override;
  };

  class TransIwR : public CloneableTransition<TransIwR> {
  public:
    TransIwR() : CloneableTransition("TransIwR") {}
    State apply(const State & ) const override;
    double rate(const State & , const Parameters & ) const override;
  };

  class TransImR : public CloneableTransition<TransImR> {
  public:
    TransImR() : CloneableTransition("TransImR") {}
    State apply(const State & ) const override;
    double rate(const State & , const Parameters & ) const override;
  };

  class TransIwH : public CloneableTransition<TransIwH> {
  public:
    TransIwH() : CloneableTransition("TransIwH") {}
    State apply(const State & ) const override;
    double rate(const State & , const Parameters & ) const override;
  };

  class TransImH : public CloneableTransition<TransImH> {
  public:
    TransImH() : CloneableTransition("TransImH") {}
    State apply(const State & ) const override;
    double rate(const State & , const Parameters & ) const override;
  };

  class TransHR : public CloneableTransition<TransHR> {
  public:
    TransHR() : CloneableTransition("TransHR") {}
    State apply(const State & ) const override;
    double rate(const State & , const Parameters & ) const override;
  };

} // namespace sarsmodel


#endif
