#include "smc.hpp"

#include <algorithm> // for_each, sort
#include <numeric> // accumulate
#include <gsl/gsl_statistics_double.h>

#include "engine.hpp"


ParticleFilter::ParticleFilter(int J, const Model & mod,
		unsigned long seed, int threads) : wp(threads) {
	rng.seed(seed);
  rel_temp = 1.0;
  effective_swarm_size = J;
  effective_sample_size = J;
	cond_log_lik_hat = 0.0;
	model = mod;
  Engine e(model.s0, model.getTransitions(), model.par);

  particles.reserve(J);
  for ( int j = 0; j < J; ++j ) {
    particles.push_back(new Particle(model.s0, model.par, e, model.obsfun, model.sde));
  }
	J_coffin = 0;
	track_coffin = false;
	// by default: don't keep track of the genealogy
	store_genealogy = false;
	// compute some initial statistics
	pred_range = getRange(ALPHA_RANGE);
	pred_median = filter_median = getMedian();
}

ParticleFilter::~ParticleFilter() {
	wp.sync();
  std::for_each(particles.begin(), particles.end(),
      [](Particle* particle){delete particle;});
}

void ParticleFilter::print(std::ostream & os) const {
	double t = 0;
	if ( !particles.empty() ) {
		t = particles.front()->getState().t();
	}
  os << "<particle_filter "
	   << "t='" << t << "' "
     << "rel_temp='" << rel_temp << "' "
     << "J_eff='" << effective_swarm_size << "' "
     << "J_inv_simpson='" << effective_sample_size << "' "
		 << "cond_LL_hat='" << cond_log_lik_hat << "' "
		 << "J_coffin='" << J_coffin << "' "
     << ">" << std::endl;
	// print statistics: prediction range:
	os << "<pred_range >" << std::endl;
	os << pred_range.first << std::endl
		 << pred_range.second << std::endl;
	os << "</pred_range>" << std::endl;
	os << "<pred_median >" << std::endl
	   << pred_median << std::endl
		 << "</pred_median>" << std::endl;
  os << "<filter_median >" << std::endl
 	   << filter_median << std::endl
		 << "</filter_median>" << std::endl;
  os << "</particle_filter>";
}

void ParticleFilter::printAllParticles(std::ostream & os) const {
  os << "<particle_filter "
     << "rel_temp='" << rel_temp << "' "
     << "J_eff='" << effective_swarm_size << "' "
     << "J_inv_simpson='" << effective_sample_size << "' "
		 << "cond_LL_hat='" << cond_log_lik_hat << "' "
		 << "J_coffin='" << J_coffin << "' "
     << ">" << std::endl;
  for ( Particle* particle : particles ) {
    os << *particle << std::endl;
  }
	// print statistics: prediction range:
	os << "<pred_range >" << std::endl;
	os << pred_range.first << std::endl
		 << pred_range.second << std::endl;
	os << "</pred_range>" << std::endl;
	os << "<pred_median >" << std::endl
	   << pred_median << std::endl
		 << "</pred_median>" << std::endl;
  os << "<filter_median >" << std::endl
 	   << filter_median << std::endl
		 << "</filter_median>" << std::endl;
  os << "</particle_filter>";
}

void ParticleFilter::reset(Rng & rng) {
  // reset all particles to state s0 (but keep their parameters)
  std::for_each(particles.begin(), particles.end(), [&](Particle* x){
		State s0j = model.initfun(x->getParameters(), rng);
		x->resetState(s0j);
	});
	J_coffin = 0;
}

void ParticleFilter::selectIndivParams(int r) {
  // reset all particles to state s0 (but keep their parameters)
  std::for_each(particles.begin(), particles.end(), [&](Particle* x){
		x->selectIndivParams(r);
	});
}

void ParticleFilter::storeGenealogy(double dt) {
	// clear pre-existing Genealogy
	G.clear();
	// set the flag
	store_genealogy = true;
	// instuct the particles to save their paths with time step dt
  std::for_each(particles.begin(), particles.end(), [&](Particle* x) {
		// remove any pre-existing path
		x->clearPath();
		// and instruct the particle to store a new path
		x->storePath(dt);
	});
}

void ParticleFilter::setCoffinTracker(bool track_coffin) {
	this->track_coffin = track_coffin;
}


std::pair<double, int> ParticleFilter::step(const Observation & obs, int J_new) {
  // if a zero or negative J_new is given, replace J_new with the current swarm size
  if ( J_new <= 0 ) { J_new = getJ(); }
  // increment paticle time and advance the state to the next time
  double t_next = obs.t; // particles will be updated to this time point
  // update the particles
  for ( Particle* particle : particles ) {
    // inform the particle of the destination time
    particle->setNextTime(t_next);
		unsigned long lseed = rng.integer();
    // and let the worker pool evolve the particle to this destination time
    wp.push_back([=](){
			// init a local RNG for the advaning the particle
			Rng lrng(lseed);
			// and use the execute method
			particle->execute(lrng);
		});
  }
  // wait for the worker pool to finish all jobs
  wp.sync();

  // not every particle will have a non-zero weight: Make a dict of those that do
  std::map<int, double> logweights;
	// keep track of the log-weights of the coffin particles
	std::map<int, double> logweights_coffin;

	// count the number of advanced coffin particles
	auto log_like_coffin = model.coflikfun(obs);

  // compute the weights of the updated particles
  std::vector<double> log_likes;
  log_likes.reserve(particles.size());
  for ( size_t j = 0; j < particles.size(); ++j ) {
		// trst if the state is the coffin state
		bool is_coffin = model.coftestfun(particles[j]->getState());
		// compute the log-likelihood
		std::pair<double, bool> loglike;
		if ( is_coffin && track_coffin ) {
			// for consistency, use the special coffin likelihood function (pre-computed)
			loglike = log_like_coffin;
		} else {
			// use the normal likelihood function
			loglike = model.likfun(obs, particles[j]->getState(), particles[j]->getParameters());
		}
    /* TODO/FIXME: the model might allow for a particular
		 * parameter hierachy: model should define a prior.
		 * TESTING: incorporate prior in the parameter random walk using
		 * a metropolis-hastings accept-reject step
		 */
    double logprior_new = particles[j]->getParameters().loglike();
    // double logprior_old = particles[j]->getCachedLogPrior(); // un-used
    // test if the particle is valid
    if ( loglike.second ) {
      double ll = loglike.first;
			// put coffin particles and normal particles in different maps
			if ( is_coffin && track_coffin ) {
    		logweights_coffin[j] = ll;
			} else {
				logweights[j] = ll;
			}
      log_likes.push_back(ll);
    } // else: do nothing
    // cache these values using the particle object
    particles[j]->cacheLogLike(loglike);
    particles[j]->cacheLogPrior(logprior_new);
  }

	// count the number of new coffin particles
	int A_coffin = logweights_coffin.size();

  // set statistics to 0 a priori
  effective_swarm_size = 0; // computed during resampling
  effective_sample_size = 0.0;

  int Jvalid = log_likes.size();
  // we cannot do re-sampling if there are no 'valid' particles
  if ( Jvalid == 0 ) {
		// filtering failure: no recovery.
  	std::cerr << "WARNING: no simulation leads to a particle with positive likelihood"
            	<< RIGHT_HERE << std::endl;
		// compute some stats...
		if ( store_genealogy ) {
			// assing to every index count one
			std::map<int, int> ones;
			for ( size_t i = 0; i < particles.size(); ++i ) { ones[i] = 1; }
			G.update(particles, ones);
		}
		pred_range = getRange(ALPHA_RANGE);
		pred_median = filter_median = getMedian();

		// apply StateModifier
		std::for_each(particles.begin(), particles.end(), [&](Particle* x) {
			model.modfun(x->getState(), obs, x->getParameters(), rng);
		});
		// after optional mods, return.
  	return std::make_pair(0.0, Jvalid); // 0.0 is a dummy value
	} // else.. resample

	if ( !track_coffin ) {
		// merge the two logweights maps
		logweights.merge(logweights_coffin);
	} else if ( logweights.empty() ) {
		// we still have coffin particles. fall back to track_coffin = false
		std::cerr << "WARNING: all particles are in the coffin state." << std::endl
		         	<< "Fall back to 'track_coffin = false' mode."
							<< RIGHT_HERE << std::endl;
		logweights.swap(logweights_coffin);
		track_coffin = false;
		J_coffin = 0; A_coffin = 0;
  }

  effective_sample_size = exp(log_eff_sam_size(log_likes));

  /* finish the estimate of the likelihood
   * also take the 'invalid' particles into account (i.e. they have ZERO likelihood)
	 * The weights of the new coffin-particles are included in log_likes.
	 * The weights of the old coffin particles (J_coffin) also
	 * have to be included.
   */
	double log_lik_tot = log_sum_exp(log_likes);
	auto & [log_w_cof, w_cof_pos] = log_like_coffin;
	if ( w_cof_pos && J_coffin > 0 ) {
		// add J_coffin times w_cof to the summed likelihood
		double log_lik_cof = log(J_coffin) + log_w_cof;
		cond_log_lik_hat = log_sum_exp(log_lik_tot, log_lik_cof) - log(getJ() + J_coffin);
	} else {
		/* coffin particles have ZERO weight: we can just use log_lik_tot,
		 * but have to include J_coffin in the denominator.
		 */
  	cond_log_lik_hat = log_lik_tot - log(getJ() + J_coffin);
	}

  // re-sample indices
  auto idxcounts = systematicResamplingLogScale(logweights, J_new, rng);

	/* Sample the number of coffin-particles.
	 * Let A be the number of non-coffin advanced particles,
	 * and A_cof the number of advanced coffin particles.
	 * We have to compute the new number of coffin particles J_cof^new.
	 * The total weight is W_tot = sum_{j=1}^A w_j + (A_cof + J_cof) w_cof
	 * which we write as W_tot = W + W_cof. We have (J + J_cof^new) * W/W_tot = J
	 * Hence, J_cof^new = J W_cof / W.
	 * We will sample J_new particles, so use J_new in the computation.
	 */
	if ( w_cof_pos && A_coffin + J_coffin > 0 ) {
		double log_W_cof = log_w_cof + log(A_coffin + J_coffin);
		double log_W = (A_coffin > 0 ? log_diff_exp(log_lik_tot, log_w_cof + log(A_coffin)) : log_lik_tot);
		double log_cof_wt_ratio = log_W_cof - log_W;
		double expect_J_coffin = J_new * exp(log_cof_wt_ratio);
		J_coffin = (expect_J_coffin > 0 ? rng.poisson(expect_J_coffin) : 0);
	} else {
		// if the coffin state has ZERO weight, they are thrown away!
		J_coffin = 0;
	}

	// keep track of the Particle Genealogy
	if ( store_genealogy ) {
		G.update(particles, idxcounts);
	}

	// before resampling, we compute some prediction statistics
	pred_range = getRange(ALPHA_RANGE);
	pred_median = getMedian();

	// replace the particles
  std::list<Particle*> resam_particles;
  for ( size_t idx = 0; idx < particles.size(); ++idx ) {
    Particle* particle = particles[idx];
    // find the index in idxcounts
    auto icit = idxcounts.find(idx);
    // when idx is not in the map, or the count is 0, the particle must be deleted
    if ( icit == idxcounts.end() || icit->second <= 0 ) {
      delete particle;
    } else {
      // add a copy of particle icit->second - 1 times
      for ( int c = 1; c < icit->second; ++c ) {
        resam_particles.push_back(new Particle(*particle));
      }
      // and add the original particle (avoid copying)
      resam_particles.push_back(particle);
      // but increase particle "generation" of the old particle
      particle->getNewId();
      effective_swarm_size += 1;
    } // if, else index count <= 0
  }
  // replace the particles with resampled particles
  particles.assign(resam_particles.begin(), resam_particles.end());

	/* Compute post-filter statistics.
	 * filter median is conditioned on not being in the "coffin" state
	 */
	filter_median = getMedian();

	/* a Model can define a StateModifier function that must be invoked
	 * in case of certain events at observation times.
	 * The default StateModifier does absolutely nothing.
	 */
	std::for_each(particles.begin(), particles.end(), [&](Particle* x) {
		model.modfun(x->getState(), obs, x->getParameters(), rng);
	});

	// return the conditional loglike and the number of valid particles
  return std::make_pair(cond_log_lik_hat, Jvalid);
}

// get functions and properties
int ParticleFilter::getJ() const { return particles.size(); }

int ParticleFilter::getEffectiveJ() const { return effective_swarm_size; }

double ParticleFilter::getEffectiveSampleSize() const { return effective_sample_size; }

int ParticleFilter::getJcoffin() const { return J_coffin; }

const Genealogy & ParticleFilter::getGenealogy() const {
	return G;
}

void ParticleFilter::coolDown(double q) {
  rel_temp *= q; // TODO: alternative cooldown schemes
  std::for_each(particles.begin(), particles.end(), [&](Particle* x){x->setRelTemp(rel_temp);});
}


std::map<int, int> systematicResampling(const std::map<int, double> & weights, Rng & rng) {
  return systematicResampling(weights, weights.size(), rng);
}

/** resampling algorithm borrowed from the pomp R package
 */
std::map<int, int> systematicResampling(const std::map<int, double> & weights, int J, Rng & rng) {
	std::map<int, int> idxcounts; // to be returned
	if ( J <= 0 ) {
    return idxcounts;
	} // else, we can assume that weights has elements
  std::list<double> sorted_weights;
	for ( auto & iw : weights ) {
		if ( iw.second < 0.0 ) {
			throw std::invalid_argument("negative weight" + RIGHT_HERE);
		} // else, add weight to the sum and add a zero to idx
		sorted_weights.push_back(iw.second);
    idxcounts[iw.first] = 0;
	}
  sorted_weights.sort();
  // now compute the cumulative weight by summing over the sorted weight list
  double sw = std::accumulate(sorted_weights.begin(), sorted_weights.end(), 0.0);
	if ( sw == 0.0 ) {
		throw std::invalid_argument("no positive weights" + RIGHT_HERE);
	}
	double dw = sw / J;
	double u = rng.uniform(0.0, dw);
	// fill idx;
  auto wit = weights.begin(); // iterator over the weights
	double cw = wit->second;
	for ( int j = 0; j < J; ++j ) {
		while ( cw < u ) {
			++wit;
			if ( wit == weights.end() ) {
				throw std::logic_error("iterator has reached end of weights" + RIGHT_HERE);
			}
			cw += wit->second;
		} // else increase the index count of the weight
		idxcounts[wit->first] += 1;
		u += dw;
	}
	return idxcounts;
}



// return the same number of indices as the number of given weights (pass log weights)
std::map<int, int> systematicResamplingLogScale(const std::map<int, double> & logweights, Rng & rng) {
  return systematicResamplingLogScale(logweights, logweights.size(), rng);
}

/** Resampling algorithm borrowed from the pomp R package
 * using weights on a log scale
 * TODO: perhaps it is just fine to use the sorted (and normalized)
 * weights on a linear scale during the final step.
 * return a given number of weights J
 */
std::map<int, int> systematicResamplingLogScale(const std::map<int, double> & logweights, int J, Rng & rng) {
	std::map<int, int> idxcounts; // to be returned
	if ( J <= 0 ) {
    return idxcounts;
	} // else, we can assume that weights has elements

  // compute the cumulative weight by summing over the sorted weight list
  std::list<std::pair<int, double>> sorted_logweights(logweights.begin(), logweights.end());
  auto comp = [](const std::pair<int, double> & p1, const std::pair<int, double> & p2) -> bool {
    return p1.second < p2.second;
  };
  sorted_logweights.sort(comp);
  //  make a vector with indices
  std::vector<int> indexvec;
  indexvec.reserve(sorted_logweights.size());
  for ( auto & iw : sorted_logweights ) {
    indexvec.push_back(iw.first);
  }

  // make a vector with cumulative weights
  std::vector<double> logcumulweights;
  logcumulweights.reserve(sorted_logweights.size());
  // add the first weight "manually"
  logcumulweights.push_back(sorted_logweights.front().second);
  for ( auto it = ++sorted_logweights.begin(); it != sorted_logweights.end(); ++it ) {
    logcumulweights.push_back(log_sum_exp(logcumulweights.back(), it->second));
  }
  // normalize the cumulative weights
  double lsw = logcumulweights.back(); // the log of the sum of the weights
  std::for_each(logcumulweights.begin(), logcumulweights.end(), [&](double & lcw){lcw -= lsw;});
  // draw a uniform(0, 1/J) distributed random number
  double u = rng.uniform(0, 1.0/J);
  size_t k = 0; // the index in the cumulative weight vector
  for ( int j = 0; j < J; ++j ) {
    double lu = log(u);
    while ( logcumulweights[k] < lu ) {
      ++k;
      if ( k >= logcumulweights.size() ) {
        throw std::logic_error("index has reached end of weights" + RIGHT_HERE);
      }
    } // else increase the index count of the weight
    idxcounts[indexvec[k]] += 1;
    u += 1.0/J;
  }
	return idxcounts;
}


std::vector<State> ParticleFilter::getQuantiles(const std::vector<double> & qs) const {
	// check validity of arguments
	for ( auto q : qs ) {
		if ( q < 0 || q > 1 ) {
			throw std::invalid_argument("range parameter must be between 0 and 1" + RIGHT_HERE);
		}
	}
	// check if the ParticleFilter object is OK
	if ( particles.empty() ) {
		throw std::logic_error("can't compute prediction range with zero particles" + RIGHT_HERE);
	}
	// construct containers for the quantiles
	State obs_zero = model.obsfun(model.s0, model.par).zero_like();
	State s_zero = model.s0.zero_like();
	// set the time (NB: particles is not empty)
	State s_obs_zero = concat(s_zero, obs_zero);
	double t = particles[0]->getState().t();
	s_obs_zero.t() = t;
	// make vector with zeros
	std::vector<State> sqs(qs.size(), s_obs_zero);
	// reserve a vector for sorting doubles (can be re-used)
	std::vector<double> xs(particles.size());
	// loop over variables
	for ( size_t i = 0; i < s_obs_zero.flat_size(); ++i ) {
		// make a vector of values
		for ( size_t j = 0; j < particles.size(); ++j ) {
			Particle* particle = particles[j];
			if ( i < s_zero.flat_size() ) {
				xs[j] = particle->getState()[i];
			} else {
				xs[j] = particle->getObservables()[i-s_zero.flat_size()];
			}
		}
		// sort the list and copy to a vector
		std::sort(xs.begin(), xs.end());
		// use GSL method to compute the quantile
		for ( size_t k = 0; k < qs.size(); ++k ) {
			sqs[k][i] = gsl_stats_quantile_from_sorted_data(xs.data(), 1, xs.size(), qs[k]);
		}
	}
	return sqs;
}

/** Compute the alpha-range of the State elements.
 */
std::pair<State, State> ParticleFilter::getRange(double alpha) const {
	// define quantiles
	double q1 = (1-alpha)/2; double q2 = 1-q1;

	std::vector<double> qs = {q1, q2};
	std::vector<State> sqs = getQuantiles(qs);
	return std::make_pair(sqs.front(), sqs.back());
}


/** Compute the median of the state */
State ParticleFilter::getMedian() const {
	// median is 0.5th quantile
	std::vector<double> q = {0.5};
	std::vector<State> sq = getQuantiles(q);
	return sq.front();
}



/** Compute the alpha-range of the State elements
 * @todo: compute quantiles in one go (incl the median
 */
std::pair<Parameters, Parameters> ParticleFilter::getParamRange(double alpha) const {
	if ( alpha < 0 || alpha > 1 ) {
		throw std::invalid_argument("range parameter must be between 0 and 1" + RIGHT_HERE);
	}
	if ( particles.empty() ) {
		throw std::logic_error("can't compute parameter range with zero particles" + RIGHT_HERE);
	}
	// define quantiles
	double f1 = (1-alpha)/2; double f2 = 1-f1;
	// construct containers for the quantiles
	Parameters par1 = model.par;
	Parameters par2 = model.par;
	// loop over parameters
	for ( size_t i = 0; i < par1.flat_size(); ++i ) {
		// only compute ranges for unlocked parameters
		if ( par1[i].isLocked() ) {
			continue;
		} // else...
		// loop over unit-specific values of the parameter (if any)
		for ( size_t r = 0; r < par1[i].size(); ++r ) {
			// make a list of values
			std::list<double> xs;
			for ( Particle* particle : particles ) {
				xs.push_back(particle->getParameters()[i][r]);
			}
			// sort the list and copy to a vector
			xs.sort();
			std::vector<double> ys(xs.size());
			ys.assign(xs.begin(), xs.end());
			// use GSL method to compute the quantile
			par1[i][r] = gsl_stats_quantile_from_sorted_data(ys.data(), 1, ys.size(), f1);
			par2[i][r] = gsl_stats_quantile_from_sorted_data(ys.data(), 1, ys.size(), f2);
		} // for loop over unit specific values
		/************ TODO: loc and scale *************/
	} // end of for loop over parameters
	return std::make_pair(std::move(par1), std::move(par2));
}


/** Compute the alpha-range of the State elements
 * @todo: compute quantiles in one go (incl the median
 */
Parameters ParticleFilter::getParamMedian() const {
	if ( particles.empty() ) {
		throw std::logic_error("can't compute parameter median with zero particles" + RIGHT_HERE);
	}
	// construct containers for the quantiles
	Parameters par = model.par;
	// loop over parameters
	for ( size_t i = 0; i < par.flat_size(); ++i ) {
		// only compute ranges for unlocked parameters
		if ( par[i].isLocked() ) {
			continue;
		} // else...
		// loop over unit-specific values of the parameter (if any)
		for ( size_t r = 0; r < par[i].size(); ++r ) {
			// make a list of values
			std::list<double> xs;
			for ( Particle* particle : particles ) {
				xs.push_back(particle->getParameters()[i][r]);
			}
			// sort the list and copy to a vector
			xs.sort();
			std::vector<double> ys(xs.size());
			ys.assign(xs.begin(), xs.end());
			// use GSL method to compute the quantile
			par[i][r] = gsl_stats_median_from_sorted_data(ys.data(), 1, ys.size());
		} // for loop over unit specific values
		/************ TODO: loc and scale *************/
	} // end of for loop over parameters
	return par;
}
