#ifndef SEIR_HPP_
#define SEIR_HPP_

#include "transition.hpp"

namespace seirmodel {

  enum VarNames {
    S=0,
    E,
    I,
    NUM_VARIABLES
  };


  enum ParSymbol {
    // index 0 for overdispersion parameter
    beta=1, // infection rate
    gamma, // E->I
    delta, // I->R
    N, // population size
    NUM_PARAMETERS
  };




  class TransSE : public CloneableTransition<TransSE> {
  public:
    TransSE() : CloneableTransition("TransSE") {};
    State apply(const State & ) const override;
    double rate(const State & , const Parameters & ) const override;
  };

  class TransEI : public CloneableTransition<TransEI> {
  public:
    TransEI() : CloneableTransition("TransEI") {};
    State apply(const State & ) const override;
    double rate(const State & , const Parameters & ) const override;
  };

  class TransIR : public CloneableTransition<TransIR> {
  public:
    TransIR() : CloneableTransition("TransIR") {};
    State apply(const State & ) const override;
    double rate(const State & , const Parameters & ) const override;
  };



} // namespace seirmodel


#endif
