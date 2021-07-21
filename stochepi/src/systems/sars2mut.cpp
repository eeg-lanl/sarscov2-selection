#include "sars2mut.hpp"

#include <algorithm> // std::max
#include "linalg.hpp" // aux for initial condition

using namespace sarsmodel;


constexpr double DAYS_IN_WEEK = 7;

// auxiliary functions

// smoothed step-function
inline double Hv(double t, double u) {
  return expit(t/u);
}

/** piecewise linear (and constant) function to model infection rates
 * with a number of breakpoints. TODO: use variable number of breakpoints
 */
double sarsmodel::infection_rate(const State & s, const Parameters & par) {
  double t = s.t();
  // uptakes
  double u1 = par[upsilon1];
  double u2 = par[upsilon2];
  double u3 = par[upsilon3];
  double u4 = par[upsilon4];
  // infection rates
  double b0 = par[beta0] * (1 - Hv(t - par[t1], u1));
  double b1 = par[beta1] * Hv(t-par[t1], u1) * (1 - Hv(t-par[t2], u2));
  double b2 = par[beta2] * Hv(t-par[t2], u2) * (1 - Hv(t-par[t3], u3));
  double b3 = par[beta3] * Hv(t-par[t3], u3) * (1 - Hv(t-par[t4], u4));
  double b4 = par[beta4] * Hv(t-par[t4], u4);
  return b0 + b1 + b2 + b3 + b4;
}

double pop_size(const State & s) {
  return s(S).value() + s(Ew).value() + s(Iw).value() + s(Em).value() + s(Im).value() + s(H).value() + s(R).value();
}


double sarsmodel::testing_rate(const State & s, const Parameters & par) {
  return par[theta];
}

// HACK: by default the import mutant freq is constant
std::function<double(double, const Parameters & par)>
    sarsmodel::import_freq_mut = [](double t, const Parameters & par){
  return par[q_mut];
};

// HACK: by default the import rate is constant
std::function<double(double, const Parameters & par)>
    sarsmodel::import_rate = [](double t, const Parameters & par){
  return par[lambda];
};



State sarsmodel::observables(const State & s, const Parameters & par) {
  State obs({"E", "I", "D", "Fm"});
  // viral load
  double Itot = s(Iw).value() + s(Im).value();
  double Etot = s(Ew).value() + s(Em).value();
  obs[0] = Etot;
  obs[1] = Itot;
  obs[2] = s(dD).value() * par[delta];
  obs[3] = (Itot > 0 ? s(Im).value() / Itot : 0.0);
  return obs;
}

State sarsmodel::gen_init_state(const Parameters & par, Rng & rng) {
  State s0({"S", "Ew", "Em", "Iw", "Im", "H", "R", "dIw", "dIm", "dD"});
  s0.t() = par[t0];
  // compute eivenvector
  std::vector<double> y0 = calc_initial_condition(par[xi], par[beta0], par[gamma],
      par[alpha], par[nu], par[omega]);
  double e0 = y0[0];
  double i0 = y0[1];
  double h0 = y0[2];

  double N = par[N0];
  double eps = par[epsilon];
  double pm = par[p_mut];
  double x = par[xi];
  // sample initial number of individuals
  s0(S) = N * (1-x) * (1-eps);
  s0(Ew) = rng.poisson(N * (1-x) * eps*(1-pm) * e0);
  s0(Em) = rng.poisson(N * (1-x) * eps*pm * e0);
  s0(Iw) = rng.poisson(N * (1-x) * eps*(1-pm) * i0);
  s0(Im) = rng.poisson(N * (1-x) * eps*pm * i0);
  s0(H) = rng.poisson(N * (1-x) * eps * h0);
  s0(R) = N*x;
  return s0;
}

// define parameter names as strings
const std::map<ParSymbol, std::string> sarsmodel::parNames = {
  std::make_pair(t0, "t0"),
  std::make_pair(N0, "N0"),
  std::make_pair(lambda, "lambda"),
  std::make_pair(q_mut, "q_mut"),
  std::make_pair(p_mut, "p_mut"),
  std::make_pair(sigma, "sigma"),
  std::make_pair(beta0, "beta0"),
  std::make_pair(beta1, "beta1"),
  std::make_pair(beta2, "beta2"),
  std::make_pair(beta3, "beta3"),
  std::make_pair(beta4, "beta4"),
  std::make_pair(t1, "t1"),
  std::make_pair(t2, "t2"),
  std::make_pair(t3, "t3"),
  std::make_pair(t4, "t4"),
  std::make_pair(upsilon1, "upsilon1"),
  std::make_pair(upsilon2, "upsilon2"),
  std::make_pair(upsilon3, "upsilon3"),
  std::make_pair(upsilon4, "upsilon4"),
  std::make_pair(alpha, "alpha"),
  std::make_pair(sarsmodel::gamma, "gamma"), // gamma is ambiguous
  std::make_pair(nu, "nu"),
  std::make_pair(omega, "omega"),
  std::make_pair(delta, "delta"),
  std::make_pair(mu, "mu"),
  std::make_pair(eta_r, "eta_r"),
  std::make_pair(delta_r, "delta_r"),
  std::make_pair(epsilon, "epsilon"),
  std::make_pair(xi, "xi"),
};



// transitions and rates

State TransSEw::apply(const State & s) const {
  State sprime = s;
  sprime(S).disc_value -= 1;
  sprime(Ew).disc_value += 1;
  return sprime;
}
double TransSEw::rate(const State & s, const Parameters & par) const {
  double betat = infection_rate(s, par);
  double N = pop_size(s);
  double q = import_freq_mut(s.t(), par);
  double lam = import_rate(s.t(), par) * std::max(0.0, 1-q); // avoid numerical issues
  double frI = (N > 0 ? s(Iw).value() / N : 0.0);
  return s(S).value() * (betat * (1-par[mu]) * frI + lam);
}

State TransSEm::apply(const State & s) const {
  State sprime = s;
  sprime(S).disc_value -= 1;
  sprime(Em).disc_value += 1;
  return sprime;
}
double TransSEm::rate(const State & s, const Parameters & par) const {
  double betat = infection_rate(s, par);
  double N = pop_size(s);
  double q = import_freq_mut(s.t(), par);
  double lam = import_rate(s.t(), par) * q;
  double frI_adj = (N > 0 ? ((1 + par[sigma]) * s(Im).value() + par[mu] * s(Iw).value()) / N : 0.0);
  return s(S).value() * (betat * frI_adj + lam);
}

State TransEwIw::apply(const State & s) const {
  State sprime = s;
  sprime(Ew).disc_value -= 1;
  sprime(Iw).disc_value += 1;
  sprime(dIw).disc_value += 1;
  return sprime;
}
double TransEwIw::rate(const State & s, const Parameters & par) const {
  return par[alpha] * s(Ew).value();
}

State TransEmIm::apply(const State & s) const {
  State sprime = s;
  sprime(Em).disc_value -= 1;
  sprime(Im).disc_value += 1;
  sprime(dIm).disc_value += 1;
  return sprime;
}
double TransEmIm::rate(const State & s, const Parameters & par) const {
  return par[alpha] * s(Em).value();
}

State TransIwR::apply(const State & s) const {
  State sprime = s;
  sprime(Iw).disc_value -= 1;
  sprime(R).disc_value += 1;
  return sprime;
}
double TransIwR::rate(const State & s, const Parameters & par) const {
  return par[gamma] * s(Iw).value();
}

State TransImR::apply(const State & s) const {
  State sprime = s;
  sprime(Im).disc_value -= 1;
  sprime(R).disc_value += 1;
  return sprime;
}
double TransImR::rate(const State & s, const Parameters & par) const {
  return par[gamma] * s(Im).value();
}


State TransIwH::apply(const State & s) const {
  State sprime = s;
  sprime(Iw).disc_value -= 1;
  sprime(H).disc_value += 1;
  return sprime;
}
double TransIwH::rate(const State & s, const Parameters & par) const {
  return par[nu] * s(Iw).value();
}

State TransImH::apply(const State & s) const {
  State sprime = s;
  sprime(Im).disc_value -= 1;
  sprime(H).disc_value += 1;
  return sprime;
}
double TransImH::rate(const State & s, const Parameters & par) const {
  return par[nu] * s(Im).value();
}

State TransHR::apply(const State & s) const {
  State sprime = s;
  sprime(H).disc_value -= 1;
  sprime(R).disc_value += 1;
  sprime(dD).disc_value += 1;
  return sprime;
}
double TransHR::rate(const State & s, const Parameters & par) const {
  return par[omega] * s(H).value();
}
