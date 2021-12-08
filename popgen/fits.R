require(rstan)
options(mc.cores=4)
rstan_options(auto_write=TRUE)
require(lubridate)

#--------------------------------------------------
# Select what to analyze
#--------------------------------------------------

# variant <- "D614G"
# datfile <- "../data/seq-pl-614_epi-17062021.csv"
# omit_m1 <- NULL
# first_day <- ymd("2020-01-04")
# cases_so_far <- NULL
# last_day <- ymd("2020-08-24")

# variant <- "B.1.351"
# datfile <- "../data/pango_2021-05-20.csv"
# omit_m1 <- "South Africa"
# first_day <- ymd("2020-10-01")
# cases_so_far <- NULL
# last_day <- ymd("2021-04-19")

# variant <- "B.1.1.7"
# datfile <- "../data/pango_2021-05-20.csv"
# omit_m1 <- "United Kingdom"
# first_day <- ymd("2020-09-20")
# cases_so_far <- NULL
# last_day <- ymd("2021-02-01")

variant <- "R.1"
datfile <- "../data/pango_2021-05-20.csv"
omit_m1 <- "Japan"
first_day <- ymd("2020-10-24")
cases_so_far <- 2500
# last_day <- ymd("2021-03-29")

remove_country <- NULL
min_num <- 20

#--------------------------------------------------
# Data for all countries.  For this variant.
#--------------------------------------------------

dat0 <- read.csv(datfile, as.is=T)

if (datfile == "../data/pango_2021-05-20.csv")
{
    dat_all <- dat0[, c("Country", "Day", "Total", variant)]
} else {
    dat_all <- dat0[, c("country", "date", "N", "N.mut")]
    dat_all <- subset(dat_all, !is.na(country))
}
names(dat_all) <- c("country", "date", "num_all", "num_new")

dat_all$freq_new <- dat_all$num_new / dat_all$num_all
dat_all$num_old <- dat_all$num_all - dat_all$num_new
dat_all$date <- ymd(dat_all$date)

if (!exists("first_day"))
    first_day <- min(subset(dat_all, num_new > 0)$date)
dat_all <- subset(dat_all, date >= first_day)
dat_all$daynum <- as.integer(dat_all$date - first_day)

#--------------------------------------------------
# Data for only the time window to analyze.
#--------------------------------------------------

if (!is.null(cases_so_far))
{
    dat_win <- subset(dat_all, num_new > 0)[, c("date", "num_new")]
    dat_win <- aggregate(num_new ~ date, data=dat_win, FUN=sum)
    dat_win <- dat_win[order(dat_win$date),]
    dat_win$cum_num_new <- cumsum(dat_win$num_new)

    last_day <- dat_win$date[min(which(dat_win$cum_num_new >= cases_so_far))]
    dat_all <- subset(dat_all, date <= last_day)
} else if (!is.null(last_day)) {
    dat_all <- subset(dat_all, date <= last_day)
    cases_so_far <- sum(dat_all$num_new)
} else {
    last_day <- NULL
    cases_so_far <- sum(dat_all$num_new)
}

#--------------------------------------------------
# Data for countries to be analyzed.
#--------------------------------------------------

if (!is.null(remove_country))
    dat_all <- subset(dat_all, !(country %in% remove_country))

places <- c()
for (place in sort(unique(dat_all$country)))
{
    dat <- subset(dat_all, country==place)
    dat <- subset(dat, num_all > 0)
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
dat$id <- sapply(dat$country, function(x) which(places == x))

#--------------------------------------------------
# Fit the hierarchical-selection-migration model.
#--------------------------------------------------

data_list <- list(N=nrow(dat), num_all=dat$num_all, num_new=dat$num_new,
                  day=dat$daynum, id=dat$id, nid=length(unique(dat$id)),
                  omit_m1=omit_m1)
stanfit <- stan(file="selection_migration.stan", data=data_list, warmup=1000,
                iter=3000, control=list(adapt_delta=0.9, max_treedepth=12))

outfile <- paste("stanfit-", variant, "-", cases_so_far, "x", min_num, ".rda", sep="")
save(stanfit, dat, places, first_day, last_day, variant, cases_so_far, min_num, file=outfile)
