/*
 * Population genetics model with selection and migration.
 * One focal variant versus all others.
 * Hierarchical across countries.
 * Transformed to logit scale.
 */

functions {

    // returns log((1-x^n)/(1-x)), dealing with the singularity at x=1
    real log_geom_sum(real logx, int n) {
        if (n == 0) {
            return negative_infinity();
        } else if (logx < -1e-5) {
            return log1m_exp(n*logx) - log1m_exp(logx);
        } else if (logx > 1e-5) {
            return log_diff_exp(n*logx, 0.0) - log_diff_exp(logx, 0.0);
        } else { // logx is very small
            return log(n) + log1p(0.5*(n-1)*expm1(logx)); // ignore H.O.T.
        }
    }
  
    // returns logit of the variant probability
    vector get_logit_prob(int N, int[] id, int[] day, vector logr0, vector logsigma, vector logmu) {
        vector[N] logit_theta;
        real logsigma_i;
        real logmu_i;
        real logr0_i;

        // the (logit) binomial probability on each day
        for (i in 1:N) {
            logsigma_i = logsigma[id[i]];
            logmu_i = logmu[id[i]];
            logr0_i = logr0[id[i]];
            
            logit_theta[i] = log_sum_exp(day[i] * logsigma_i + logr0_i,
                logmu_i + log_geom_sum(logsigma_i, day[i])
            );
        }

        return logit_theta;
    }
}

data {
    int<lower=1> N;                // number of data points (all countries combined)
    int<lower=0> day[N];           // day number (choose the 0 point to suit the variant considered)
    int<lower=1> num_all[N];       // number of 'trials'/samples/sequences each day
    int<lower=0> num_new[N];       // number of 'successes'/new variant samples each day
    int<lower=1> nid;              // number of countries
    int<lower=1, upper=nid> id[N]; // country corresponding to each data point
    int<lower=0, upper=1> omit_m1; // whether to omit migration into the first country
}

parameters {
    vector[nid] logsigma;                    // log((1+s)/(1-m))
    vector[nid] logr0;                       // logit(p_0)
    vector[omit_m1 ? nid-1 : nid] logmu_raw; // log(m/(m-1)) = logit(m)
  
    // for normal prior on s 
    real s_mean;
    real<lower=0> s_sd;
}

transformed parameters {
    vector<lower=-1>[nid] s;
    vector<lower=0, upper=1>[nid] p0;
    vector[nid] logmu;
    vector<upper=1>[nid] m; // (should have lower=0 when omit_m1=FALSE)

    // omit migration into the first country, if specified
    if (omit_m1) {
        logmu = append_row(negative_infinity(), logmu_raw);
    } else {
        logmu = logmu_raw;
    }
    m = inv_logit(logmu);
    s = exp(logsigma) .* (1.0 - m) - rep_vector(1.0, nid);
    p0 = inv_logit(logr0);
}

model {

    /*** model ***/
 
    // get the binomial probability on each day
    vector[N] logit_theta = get_logit_prob(N, id, day, logr0, logsigma, logmu);
    // put into the likelihood
    num_new ~ binomial_logit(num_all, logit_theta);

    /*** priors ***/

    // separate initial frequency for each country
    p0 ~ uniform(0, 1);

    // migration from a specified distribution
    if (omit_m1) {
        for (i in 2:nid) {
            m[i] ~ exponential(1000);
        }
    } else {
        m ~ exponential(1000);
    }

    // selection from a normal distribution, estimated in the hierarchical model
    s ~ normal(s_mean, s_sd);
    s_mean ~ normal(0, 0.5);      // centered on 0 for no selection
    s_sd ~ normal(0.25, 0.25);    // just seems reasonable

    /*** jacobian corrections ***/

    target += 2 * log(p0) - logr0;          // log|dp0/dlogr0|
    target += logsigma + log1m(m);          // log|ds/dlogsigma|
    if (omit_m1) {
        for (i in 2:nid) {
            target += 2 * log(m[i]) - logmu[i];
        }
    } else {
        target += 2 * log(m) - logmu;       // log|dm/dlogmu|
    }
}

generated quantities {
    vector[nid] mgen;
    vector[nid] sgen;
    real sgen_mean;
    real sgen_sd;
    real gentime = -1;
    vector[N] theta;

    /*** convert to units of per generation ***/

    while (gentime < 0) {
        // chosen so mean serial interval = 4-7.8 days is plausible (and not negative)
        gentime = normal_rng(5.9, 1.15);
    }

    for (i in 1:nid) {
        mgen[i] = pow(1 + m[i], gentime) - 1;
        sgen[i] = pow(1 + s[i], gentime) - 1;
    }
    sgen_mean = pow(1 + s_mean, gentime) - 1;
    sgen_sd   = pow(1 + s_sd,   gentime) - 1;

    /*** predicted probabilities ***/

    theta = inv_logit(get_logit_prob(N, id, day, logr0, logsigma, logmu));
}

