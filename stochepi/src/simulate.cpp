#include "simulate.hpp"

#include <fstream> // disk IO
#include <sstream>
#include <iostream>
#include <string>

#include "parallelism.hpp"
#include "rng.hpp"
#include "engine.hpp"


class ModelJob {
public:
  ModelJob(const Model & model, double tmax, double dt) :
      model(model), tmax(tmax), dt(dt) { /* empty */ }
  // result to write to file when all Jobs are finished
  std::string result;
  void execute(Rng & rng) {
    Engine engine(model.s0, model.getTransitions(), model.par);
    std::stringstream ss;
    ss << "<sim model_name='" << model.getName() << "' >" << std::endl;
    // initial time and state
    State s = model.initfun(model.par, rng);
    State obs = model.obsfun(s, model.par);
    ss << concat(s, obs) << std::endl;
    while ( s.t() < tmax ) {
      double tnext = std::min(tmax, s.t() + dt);
      s = engine.evolve_hybrid(s, tnext, model.par, rng, model.sde);
      State obs = model.obsfun(s, model.par);
      ss << concat(s, obs) << std::endl;
    }
    ss << "</sim>" << std::endl;
    result = ss.str();
  }
protected:
  const Model & model;
  double tmax;
  double dt;
};


void simulate(const Model & model, double tmax, double dt, size_t J,
      std::ostream & stream, int threads, unsigned long seed) {
  // init Rng and WorkerPool
  Rng rng(seed);
  WorkerPool wp(threads);

  // Use ModelJob to store the output of simulations
  std::list<ModelJob*> jobs;
  for ( size_t j = 0; j < J; ++j ) {
    // Use a ModelJob object.
    ModelJob* job = new ModelJob(model, tmax, dt);
    jobs.push_back(job);
    // sample local seed with global RNG
    unsigned long lseed = rng.integer();
    wp.push_back([=](){
      // use local seed to seed a local RNG
      Rng lrng(lseed);
      // and execute the job using the local RNG
      job->execute(lrng);
    });
  }
  // sync the WorkerPool before writing results to the stream.
  wp.sync_timed();

  // write results to file
  stream << "<sims model_name='" << model.getName() << "'>" << std::endl;
  for ( auto job : jobs ) {
    stream << job->result;
  }
  stream << "</sims>" << std::endl;

  // clean-up the mess
  std::for_each(jobs.begin(), jobs.end(), [](ModelJob* job){delete job;});
}



void simulate(const Model & model, double tmax, double dt, size_t J,
      std::string filename, int threads, unsigned long seed) {
  // write results to file
  std::ofstream file;
  file.open(filename.c_str());
  if ( !file.is_open() ) {
    throw std::runtime_error("could not open output file" + RIGHT_HERE);
  }
  // use the more general simulate function
  simulate(model, tmax, dt, J, file, threads, seed);
  // close the file before returning
  file.close();
}
