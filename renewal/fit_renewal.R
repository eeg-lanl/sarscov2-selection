require(tidyverse)
require(lubridate)
require(cmdstanr)
require(loo)
require(countrycode)

#--------------------------------------------------
# Data
#--------------------------------------------------

# The choice of data was all made in data_renewal_omicron.R.
# Now format it for Stan.

filedate <- "mar11"
variant <- "omicron"

# from data_renewal.R
daily_data <- read_csv(sprintf("data/%s-day-%s.csv", variant, filedate)) %>%
              arrange(Country)

# Stan won't run with this country.  Probably because some AggrCases are 0.
if (variant == "delta")
    daily_data <- daily_data %>% filter(Country != "Papua New Guinea")

filedate <- "mar11H"
i_cvr <- NA
# for covariates, see section below
# filedate <- "mar11G"
# i_cvr <- 5

# number of days in the time period for each region
num_obs <- daily_data %>% group_by(Country) %>% summarize(m = max(daynum)+1) %>% deframe()
regions <- names(num_obs)

# first day for each region
start_dates <- daily_data %>% select(Country, first_day) %>% unique() %>% deframe()
stopifnot(all(names(start_dates) == regions))

# average cases per day for each country
mean_cases <- daily_data %>%
              group_by(Country) %>%
              summarize(scaling = mean(Cases, na.rm=T))
# (if any days happen to be missing, would be better to use the total time
# duration as the denominator, under the assumption that missed days of
# cases are reported another day)

# matrix of (scaled) daily case counts
DayCases <- daily_data %>%
            select(Country, Cases, daynum) %>%
            left_join(mean_cases) %>%
            mutate(ScaledCases = Cases / scaling) %>%
            select(-Cases, -scaling) %>%
            pivot_wider(names_from=Country, values_from=ScaledCases) %>%
            select(-daynum)

# matrices of daily total sequences and variant sequences
DayTotalSeq <- daily_data %>%
               select(Country, N_total, daynum) %>%
               pivot_wider(names_from=Country, values_from=N_total) %>%
               select(-daynum)
DayVariantSeq <- daily_data %>%
               select(Country, N_focal, daynum) %>%
               pivot_wider(names_from=Country, values_from=N_focal) %>%
               select(-daynum)

# Don't actually need DayCases, DayTotalSeq, or DayVariantSeq.  Their aggregated
# versions, below, are what's used by Stan.  So could just start with
# variant-week-filedate.csv, divide weekly mean cases by 7, and set num_obs
# NumObs ObsTimes from week info.  But DayCases is used to set an initial
# condition, and they're used for platting.

#--------------------------------------------------
# Aggregate, from day to week
#--------------------------------------------------

# days per week
obs_interval <- 7

# weeks per country
NumObs <- num_obs[regions] %/% obs_interval

# day numbers corresponding to weeks
ObsTimes <- sapply(regions, function(x) seq(obs_interval, num_obs[x], obs_interval))

# aggregate everything by week
weekly_data <- daily_data %>%
               mutate(weeknum = daynum %/% 7) %>%
               group_by(weeknum, Country, first_day) %>%
               summarize(N_focal = sum(N_focal), N_other = sum(N_other),
                         N_total = sum(N_total), Cases = sum(Cases, na.rm=T)) %>%
               ungroup() %>%
               mutate(Freq_focal = N_focal / N_total)

# matrix of (scaled) weekly case counts
AggrCases <- weekly_data %>%
             select(Country, Cases, weeknum) %>%
             left_join(mean_cases) %>%
             mutate(ScaledCases = Cases / scaling) %>%
             select(-Cases, -scaling) %>%
             pivot_wider(names_from=Country, values_from=ScaledCases) %>%
             select(-weeknum)

# matrices of weekly total sequences and variant sequences
AggrTotalSeq <- weekly_data %>%
                select(Country, N_total, weeknum) %>%
                pivot_wider(names_from=Country, values_from=N_total) %>%
                select(-weeknum)
AggrVariantSeq <- weekly_data %>%
                select(Country, N_focal, weeknum) %>%
                pivot_wider(names_from=Country, values_from=N_focal) %>%
                select(-weeknum)

#--------------------------------------------------
# Covariates
#--------------------------------------------------

# C = ED10
# D = ED20
# E = VaxFull
# F = VaxBoost
# G = frac_delta

# excess deaths, as a proxy for natural infection/immmunity
cvr_ed <- read_csv("data/excess_deaths_per100k_perday_2021_04_01.csv") %>%
          select(Country = country, ED10 = ed_per100k_perday_10, ED20 = ed_per100k_perday_20)

# vaccination
# for each country, use the most recent day (up to Dec 15) with data
cvr_vx <- read_delim("data/vac_01-06-22.csv", delim=" ") %>% select(-i) %>%
          select(ISO = location, Day = date, VaxFull = percent_fullvac, VaxBoost = percent_boosted) %>%
          mutate(Country = countrycode(ISO, origin="iso3c", destination="country.name")) %>%
          filter(Day <= ymd("2021-12-15")) %>%
          drop_na() %>%
          group_by(Country) %>%
          filter(Day == max(Day)) %>%
          ungroup() %>%
          filter(Country %in% regions) %>%
          select(-ISO, -Day)

# amount of delta
# the proportion of Delta in each country, for the first week before Omicron arrived
cvr_dl <- read_csv("data/delta-2022-03-11.csv") %>%
          mutate(Country=countrycode(Country, origin="country.name", destination="country.name")) %>%
          filter(Country %in% regions) %>%
          left_join(unique(select(daily_data, Country, first_day))) %>%
          filter(Day < first_day & Day >= first_day-7) %>%
          group_by(Country) %>%
          summarize(n_delta = sum(Delta), n_other = sum(Other), frac_delta = n_delta / (n_delta + n_other)) %>%
          select(-n_delta, -n_other)

cvr_all <- inner_join(cvr_ed, cvr_vx) %>% inner_join(cvr_dl) %>%
           arrange(match(Country, regions))
stopifnot(all(cvr_all$Country == regions))

# pick one, and standardize for stan
if (!is.na(i_cvr))
    cvr1 <- scale(structure(cvr_all[[i_cvr+1]], names = cvr_all$Country))

#--------------------------------------------------
# Fit the Stan model, and save the results
#--------------------------------------------------

R <- length(regions)

if (variant == "omicron") {
    mu_kernel_wt = 4.0
    mu_kernel_mt = 2.0
} else if (variant == "delta") {
    mu_kernel_wt = 6.0
    mu_kernel_mt = 4.0
}

data_list <- list(
  R = R,
  N = num_obs,
  NumObs = NumObs,
  Cases = as.matrix(AggrCases),
  TotalSeq = as.matrix(AggrTotalSeq),
  VariantSeq = as.matrix(AggrVariantSeq),
  ObsTimes = ObsTimes,
  mu_kernel_wt = mu_kernel_wt,
  mu_kernel_mt = mu_kernel_mt,
  phi_kernel = 4,
  gentime_distr = 1,
  W = 15,
  NumCovariates = 0,
  Covariates = matrix(0.0, nrow=R, ncol=0)
)
if (!is.na(i_cvr))
{
    data_list[[NumCovariates]] <- 1
    data_list[[Covariates]] <- matrix(cvr1, nrow=R, ncol=1)
}

init_fun <- function() {
  list(
    Z_concat = array(0, sum(num_obs[regions]) - 2*R),
    tau_pop = 0.01,
    Rt0 = rep(1.0, R),
    dI0 = colMeans(DayCases[1:3,], na.rm=T) + 0.1,
    b0 = rep(1e-2, R),
    fitness = rep(1.5, R),
    mu_fitness = 0.0,
    sd_fitness = 1.0,
    covariate_weights = array(0.0, dim=1)
  )
}
# some initial parameters are still missing: consider adding them?

sm <- cmdstan_model("focalvariant-renewal-aggr.stan")

# run stan
sam <- sm$sample(
  data = data_list, 
  init = init_fun, 
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000, 
  adapt_delta = 0.95,
  max_treedepth = 13,
  #fixed_param = TRUE,
  output_dir = "stan-cache",
  refresh = 10
)

# use rstan to get MCMC output in nice format
sam.rstan <- rstan::read_stan_csv(sam$output_files())
chain <- rstan::extract(sam.rstan)

# compute LOO-IC
log_lik <- extract_log_lik(sam.rstan, parameter_name = "loglik")
loo <- loo(log_lik)
print(loo)

# write data and fit to results/
save(weekly_data, daily_data, DayCases, DayTotalSeq, DayVariantSeq, data_list, chain, regions, variant, mean_cases, start_dates, cvr_all, file = paste0("results/", variant, "-", filedate, "-ren-stanout.rda"))
