#include "seir.hpp"

using namespace seirmodel;


State TransSE::apply(const State & s) const {
  State sprime = s;
  sprime(S).disc_value -= 1;
  sprime(E).disc_value += 1;
  return sprime;
}
double TransSE::rate(const State & s, const Parameters & par) const {
  return par[beta] * s(S).value() * s(I).value() / par[N];
}

State TransEI::apply(const State & s) const {
  State sprime = s;
  sprime(E).disc_value -= 1;
  sprime(I).disc_value += 1;
  return sprime;
}
double TransEI::rate(const State & s, const Parameters & par) const {
  return par[gamma] * s(E).value();
}

State TransIR::apply(const State & s) const {
  State sprime = s;
  sprime(I).disc_value -= 1;
  // don't model R explicitly
  return sprime;
}
double TransIR::rate(const State & s, const Parameters & par) const {
  return par[delta] * s(I).value();
}
