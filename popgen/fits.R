require(rstan)
options(mc.cores=parallel::detectCores()-2)
rstan_options(auto_write=TRUE)
require(lubridate)

#--------------------------------------------------
# Select which variant to analyze
#--------------------------------------------------

### B.1.1.7 ###

datfile <- "../../data/seq-pl_epi-26022021.csv"
num_all <- "N"
num_new <- "N_mt"
min_num <- 5
omit_m1 <- "United Kingdom"
outfile <- "fits-b117.rda"

### D614G ###

datfile <- "../../data/seq_epi-08282020.csv"
num_all <- "seq_total"
num_new <- "seq_G"
min_num <- 10
omit_m1 <- NULL
outfile <- "fits-614.rda"

#--------------------------------------------------
# Data for all countries
#--------------------------------------------------

dat0 <- read.csv(datfile, as.is=T)
dat_all <- data.frame(country=dat0$country, date=dat0$date, num_all=dat0[,num_all], num_new=dat0[,num_new])
rm(dat0)

dat_all$freq_new <- dat_all$num_new / dat_all$num_all
dat_all$num_old <- dat_all$num_all - dat_all$num_new
dat_all$date <- ymd(dat_all$date)

dat_all <- subset(dat_all, !is.na(dat_all$num_all))
dat_all <- subset(dat_all, !is.na(dat_all$country))
dat_all <- subset(dat_all, date >= ymd("2019-10-01"))

#--------------------------------------------------
# Data for countries to be analyzed
#--------------------------------------------------

places <- c()
for (place in sort(unique(dat_all$country)))
{
    dat <- subset(dat_all, country==place)      # this country
    dat <- subset(dat, num_all > 0)             # days with any samples
    # "cov.lanl.gov uses minimum of 10 sequences representing a variant in the virus, with at least 14 days of sampling"
    if (nrow(dat) < 14 | sum(dat$num_old) < min_num | sum(dat$num_new) < min_num)
        next
    places <- c(places, place)
}

if (!is.null(omit_m1))
{
    places <- c(omit_m1, places[-which(places == omit_m1)])
    omit_m1 <- TRUE
} else {
    omit_m1 <- FALSE
}

dat <- subset(dat_all, country %in% places)
dat <- subset(dat, num_old > 0 | num_new > 0)

first_day <- min(subset(dat, num_new > 0)$date)
dat <- subset(dat, date >= first_day)
dat$daynum <- as.integer(dat$date - first_day)

dat$id <- sapply(dat$country, function(x) which(places == x))

data_list <- list(N=nrow(dat), num_all=dat$num_all, num_new=dat$num_new, day=dat$daynum, id=dat$id, nid=length(unique(dat$id)), omit_m1=omit_m1)

#--------------------------------------------------
# Fit the hierarchical model
#--------------------------------------------------

fitN <- stan(file="selection_migration.stan", data=data_list, warmup=1000, iter=2000, control=list(adapt_delta=0.9, max_treedepth=12))

save(fitN, dat_all, dat, places, first_day, file=outfile)
