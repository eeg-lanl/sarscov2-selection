data {
    int<lower=0> N;
    int<lower=0> NumSam[N];
    int<lower=0> NumVar[N];
    vector[N] T;
    real T_G;
}
parameters {
    real alpha;
    real beta;
}
model {
    vector[N] logit_p = alpha * T + beta;
    NumVar ~ binomial_logit(NumSam, logit_p);
}
generated quantities {
    //real s = exp(alpha*T_G) - 1;
    real s = alpha*T_G;
    vector[N] phat = inv_logit(alpha * T + beta);
}
