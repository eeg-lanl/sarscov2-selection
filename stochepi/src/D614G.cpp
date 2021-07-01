#include "D614G.hpp"

#include <fstream> // disk IO
#include <sstream> // save intermediate output
#include <algorithm> // for_each, ...

#include "rng.hpp"
#include "variable.hpp"
#include "state.hpp"
#include "macros.hpp"
#include "systems/sars2mut.hpp"
#include "model.hpp"
#include "data.hpp"
#include "parse_params.hpp"
#include "ipf.hpp"
#include "simulate.hpp"
#include "aux.hpp" // eigenvectors of 2x2 matrices
#include "splines.hpp"


void simulate_sars2mut_model(unsigned long seed, int threads, int J,
    const std::string & paramFileName, double tmax) {
  Rng rng(seed);

  // import paramter values
  std::list<ParSpecs> parSpecsList = parse_param_file(paramFileName);
  Parameters par(sarsmodel::NUM_PARAMETERS);
  load_param_specs<sarsmodel::ParSymbol>(par, sarsmodel::parNames, parSpecsList);

  // import splines
  std::string linfunFileNameFreq = "Fm_lin_fun.txt"; // FIXME!
  splines::PiecewiseLinFunc plfun_freq = splines::import_piecewise_lin_fun(linfunFileNameFreq);
  // replace default []()
  sarsmodel::import_freq_mut = [plfun_freq](double t, const Parameters & par) -> double {
    return plfun_freq(t);
  };

  std::string linfunFileNameRate = "I_lin_fun.txt"; // FIXME!
  splines::PiecewiseLinFunc plfun_rate = splines::import_piecewise_lin_fun(linfunFileNameRate);
  // replace default []()
  sarsmodel::import_rate = [plfun_rate](double t, const Parameters & par) -> double {
    // set closing time of borders
    double t_close = 355; // FIXME!
    // scale given rates spline with parameter lambda
    return ( t < t_close ? plfun_rate(t) * par[sarsmodel::lambda] : 0.0);
  };


  // overdispersion of intrinsic noise
  par[IDX_SIGMA_OD] = 0.02;

  // "S", "Ew", "Em", "Iw", "Im", "H", "R", "dIw", "dIm", "dD"
  par.select_od = {1, 1, 1, 1, 1, 1, 1, 0, 0, 0};

  // mutate the special overdispersion parameter
  par[IDX_SIGMA_OD].setNameBoundsPstdAndUnlock("sigma_od", 0.0, 0.1, 1e-4);

  double dt = 1.0;

  // define model transitions
  std::vector<Transition*> transitions = {
    new sarsmodel::TransSEw,
    new sarsmodel::TransSEm,
    new sarsmodel::TransEwIw,
    new sarsmodel::TransEmIm,
    new sarsmodel::TransIwR,
    new sarsmodel::TransImR,
    new sarsmodel::TransIwH,
    new sarsmodel::TransImH,
    new sarsmodel::TransHR,
  };

  State s0 = sarsmodel::gen_init_state(par, rng);

  // create a model object
  Model model(s0, transitions, par, true, "sars_model");

  // sample the initial number of infected cells from a poisson distribution
  model.initfun = &sarsmodel::gen_init_state;
  // additional output of the model
  model.obsfun = &sarsmodel::observables;

  std::string filename = "data/out/sars-model-sims.xml";

  std::cout << "starting simulations..." << std::endl;

  simulate(model, tmax, dt, J, filename, threads, seed);

  // clean-up the mess
  std::for_each(transitions.begin(), transitions.end(),
      [](Transition* trans){delete trans;});
}



void filter_sars2mut_model(unsigned long seed, int threads, int J, int G, int M,
    int D, const std::map<std::string, double> & fixedParamMap,
    const std::string & dataFileName, const std::string & paramFileName,
    const std::string & id) {
  // call more general method with empty filenames for the splines
  filter_sars2mut_model(seed, threads, J, G, M, D, fixedParamMap,
      dataFileName, paramFileName, id, false, "", false, 0.0);
}


// general method called by simplified functions
void filter_sars2mut_model(unsigned long seed, int threads, int J, int G, int M,
    int D, const std::map<std::string, double> & fixedParamMap,
    const std::string & dataFileName, const std::string & paramFileName,
    const std::string & id, bool timeDepMigr, const std::string & migrFileName,
    bool endTimeProvided, double endTime) {
  Rng rng(seed);

  // import data
  std::vector<Timeseries> txss = load_vector_panel_data(dataFileName);
  // if an end time is provided, prune the timeseries
  if ( endTimeProvided ) {
    for ( auto & txs : txss ) {
      pruneTimeseries(txs, endTime);
    }
  }

  for ( auto & txs : txss ) {
    std::cout << txs << std::endl; // take a look at the panel
  }

  // data starts on a Tuesday, so start sims on a Monday to complete the week!
  int R = txss.size(); // number of 'replicates'

  // import data from given files

  if ( timeDepMigr ) {
    // try loading the piecewise linear functions from file
    std::list<splines::PiecewiseLinFunc> funlist = splines::import_multiple_piecewise_lin_funs(migrFileName);

    if ( funlist.size() == 2 ) {
      splines::PiecewiseLinFunc plfun_freq = funlist.front();
      // replace default []()
      sarsmodel::import_freq_mut = [plfun_freq](double t, const Parameters & par) -> double {
        return plfun_freq(t);
      };

      splines::PiecewiseLinFunc plfun_rate = funlist.back();
      // replace default []()
      sarsmodel::import_rate = [plfun_rate](double t, const Parameters & par) -> double {
        // set closing time of borders
        double t_close = 355; // HARDCODED CLOSING DATE: FIXME!
        // scale given rates spline with parameter lambda
        return ( t < t_close ? plfun_rate(t) * par[sarsmodel::lambda] : 0.0);
      };
    } else {
      throw std::runtime_error("unable to import external FOI functions" + RIGHT_HERE);
    }
  }

  //////////////////////////////////////////////////////////////////////

  // import parameter values
  std::list<ParSpecs> parSpecsList = parse_param_file(paramFileName);
  Parameters par(sarsmodel::NUM_PARAMETERS);
  load_param_specs<sarsmodel::ParSymbol>(par, sarsmodel::parNames, parSpecsList);

  // overdispersion of intrinsic noise

  // "S", "Ew", "Em", "Iw", "Im", "H", "R", "dIw", "dIm", "dD"
  par.select_od = {1, 1, 1, 1, 1, 1, 1, 0, 0, 0};

  par[IDX_SIGMA_OD] = 0.02;
  // mutate the special overdispersion parameter
  par[IDX_SIGMA_OD].setNameBoundsPstdAndUnlock("sigma_od", 0.0, 0.05, 1e-4);

  // some parameters should only random walk during parts of the timeseries
  double u_beta = 7;
  double u_t = 30;
  double t0 = par[sarsmodel::t0];
  double t1 = par[sarsmodel::t1];
  double t2 = par[sarsmodel::t2];
  double t3 = par[sarsmodel::t3];
  // only mutate when infection rate is relavant
  par[sarsmodel::beta0].setRWRule([=](double t){return t < t1+u_beta;});
  par[sarsmodel::beta1].setRWRule([=](double t){return t > t1-u_beta && t < t2+u_beta;});
  par[sarsmodel::beta2].setRWRule([=](double t){return t > t2-u_beta && t < t3+u_beta;});
  par[sarsmodel::beta3].setRWRule([=](double t){return t > t3-u_beta;});
  // only mutate around the break points
  par[sarsmodel::t1].setRWRule([=](double t){return t > t1 - u_t && t < t1 + u_t;});
  par[sarsmodel::t2].setRWRule([=](double t){return t > t2 - u_t && t < t2 + u_t;});
  par[sarsmodel::t3].setRWRule([=](double t){return t > t3 - u_t && t < t3 + u_t;});
  // same for uptake
  par[sarsmodel::upsilon1].setRWRule([=](double t){return t > t1 - u_t && t < t1 + u_t;});
  par[sarsmodel::upsilon2].setRWRule([=](double t){return t > t2 - u_t && t < t2 + u_t;});
  par[sarsmodel::upsilon3].setRWRule([=](double t){return t > t3 - u_t && t < t3 + u_t;});
  // only mutate at time t0
  par[sarsmodel::epsilon].setRWRule([=](double t){return t <= t0;});
  par[sarsmodel::p_mut].setRWRule([=](double t){return t <= t0;});

  // for computing a likelihood profile, we want to lock some of the parameters
  std::stringstream file_id_stream;
  if ( !id.empty() ) {
    file_id_stream << "_" << id;
  }
  for ( auto & [parname, parval] : fixedParamMap ) {
    par[parname] = parval;
    par[parname].lock();
    std::cout << "fixed " << parname << " to " << parval << std::endl;
    file_id_stream << "_" << parname << "=" << parval;
  }
  std::string file_id = file_id_stream.str(); // added to model name

  // define the model
  std::vector<Transition*> transitions = {
    new sarsmodel::TransSEw,
    new sarsmodel::TransSEm,
    new sarsmodel::TransEwIw,
    new sarsmodel::TransEmIm,
    new sarsmodel::TransIwR,
    new sarsmodel::TransImR,
    new sarsmodel::TransIwH,
    new sarsmodel::TransImH,
    new sarsmodel::TransHR,
  };

  // initial state
  State s0 = sarsmodel::gen_init_state(par, rng);

  // create a model object
  Model model(s0, transitions, par, true, "sars_model" + file_id);

  // sample the initial number of infected cells from a poisson distribution
  model.initfun = &sarsmodel::gen_init_state;

  // additional output of the model
  model.obsfun = &sarsmodel::observables;

  // likelihood function
  model.likfun = [=](const Observation & obs, const State & s,
      const Parameters & par) -> std::pair<double, bool> {
    double ll = 0.0; bool ok = true;
    // likelihood of sequence samples
    double dItot = s(sarsmodel::dIw).value() + s(sarsmodel::dIm).value();
    int seq_total = obs.x[2];
    if ( dItot > 0 ) {
      double Fm = s(sarsmodel::dIm).value() / dItot;
      int seq_mut = obs.x[1];
      double r = par[sarsmodel::eta_r];
      auto [ll_seq, ok_seq] = beta_binomial_lpmf(seq_mut, seq_total, Fm, r); // k, n, p, r
      if ( ok_seq ) {
        ll += ll_seq;
      } else {
        ok = false;
      }
    } else if ( seq_total > 0 ) { // no infectious individuals means no samples possible
      ok = false;
    } // else ok = true, ll = 0.0
    // likelihood of deaths
    int deaths_cc = obs.c[0];
    if ( deaths_cc == UNCENSORED_CODE ) {
      int deaths = obs.x[0];
      double dD = std::max(0.0, s(sarsmodel::dD).value());
      double r = par[sarsmodel::delta_r];
      double delta = par[sarsmodel::delta];
      auto [ll_deaths, ok_deaths] = neg_binomial_lpmf(deaths, dD*delta, r);
      //auto [ll_deaths, ok_deaths] = poisson_lpmf(deaths, dD*delta);
      if ( ok_deaths ) {
        ll += ll_deaths;
      } else {
        ok = false;
      }
    }

    // return the total ll
    return std::make_pair(ll, ok);
  };

  // after each case observation, we have to reset the cumulative incidence
  model.modfun = [](State & s, const Observation & obs, const Parameters & par,
      Rng & rng) {
    if ( obs.event == "[RESET_CASES]" ) {
      s(sarsmodel::dIw) = int(0); // reset infectiousness incidence
      s(sarsmodel::dIm) = int(0); // idem for mutant
      s(sarsmodel::dD) = int(0); // reset confirmed deaths
    }
  };

  // clean-up the mess
  std::for_each(transitions.begin(), transitions.end(),
      [](Transition* trans){delete trans;});

  // now we are ready to start the particle filter algorithm (defined in ipf.cpp)
  AlgoParameters algpar = {J, M, D, G, 0.25};
  iterated_filtering(model, txss, algpar, threads, seed);
}
