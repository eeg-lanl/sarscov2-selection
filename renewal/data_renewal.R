# Prepare data for the renewal model: sequence counts, case counts
# Decisions: which countries, which dates for each country

require(tidyverse)
require(lubridate)
require(countrycode)
require(logihist)
require(patchwork)

#--------------------------------------------------
# Prepare variant data
#--------------------------------------------------

variant <- "omicron"
filedate <- "mar11"

# Sequence counts
dat0 <- read_csv(file.path("data", paste0(variant, "-2022-03-11.csv"))) %>%
        mutate(Day = ymd(Day)) %>%
        rename(N_focal = str_to_title(variant), N_other = Other) %>%
        mutate(N_total = N_focal + N_other, Freq_focal = N_focal / N_total)

# Clean up country names
dat0 <- dat0 %>%
        mutate(Country = countrycode(
            countrycode(Country, origin="country.name", destination="iso3c"),
            origin="iso3c", destination="country.name"))

# Some "Country" are not really countries
# (look out for more of these in places() below)
# (Aruba is a country, but we have no case count data from there.)
noncountry <- c("Curaçao", "French Guiana", "Gibraltar", "Hong Kong SAR China",
                "Northern Mariana Islands", "Puerto Rico", "Réunion", "U.S.
                Virgin Islands", "Aruba")
dat0 <- dat0 %>%
        filter(!(Country %in% noncountry))

# Criteria for "enough" data for a country
min_day <-   21  # at least this many days with data (will be further adjusted below)
min_seq <- 1000  # at least this many counts of the variant (and other)
places <- c()
for (place in sort(unique(dat0$Country)))
{
    dat1 <- dat0 %>% filter(Country == place)

    # skip if not enough data
    if (nrow(dat1) < min_day |
        sum(dat1$N_other) < min_seq |
        sum(dat1$N_focal) < min_seq)
        next

    # otherwise, include
    places <- c(places, place)
}

# Data for just these countries
dat1 <- dat0 %>%
        filter(Country %in% places)

# check: counts per country
dat1 %>% group_by(Country) %>%
         summarize(focal = sum(N_focal), other = sum(N_other)) %>%
         arrange(focal) %>% data.frame

# First day used is when X cases and Y% of cases have accumulated for the variant.
# And variant freq is at least Z%.
# (goal is to skip the stochastic left tail)
first_day_each <- dat1 %>%
                  group_by(Country) %>%
                  mutate(cs = cumsum(N_focal)) %>%
                  mutate(csf = cs / sum(N_focal)) %>%
                  filter(cs >= 10 & csf >= 0.0005 & Freq_focal > 0.1) %>%
                  summarize(first_day = min(Day))
dat2 <- left_join(dat1, first_day_each) %>%
        filter(Day >= first_day)

#--------------------------------------------------
# Prepare case count data
#--------------------------------------------------

# Being sure to keep days that have cases but not seqs, and vice versa

cas0 <- read_csv("data/epi-data.csv")
# all(places %in% cas0$country) # TRUE
cas1 <- cas0 %>%
        mutate(Day = ymd(date)) %>%
        select(Country = country, Cases = cases, Day) %>%
        mutate(Cases = replace(Cases, which(Cases < 0), NA))
cas2 <- cas1 %>%
        filter(Country %in% places) %>%
        left_join(first_day_each) %>%
        filter(Day >= first_day)

dat3 <- full_join(dat2, cas2)
#filter(dat3, is.na(Cases)) # only a few, and recent
#filter(dat3, is.na(N_total)) # several hundred

dat3 <- dat3 %>%
        replace_na(list(N_focal=0, N_other=0, N_total=0))

# Number the days
dat3 <- dat3 %>%
        mutate(daynum = as.integer(Day - first_day))

#--------------------------------------------------
# Plots
#   variant logihist
#   case counts
#   goal: sanity check, and decide on last_day
#--------------------------------------------------

# (for logihist)
convert_to_bernoulli <- function(x)
{
    val <- c("old"=0, "new"=1)
    variant <- c(rep(val["old"], x[1]), rep(val["new"], x[2]))
    if (length(variant) > 0)
    {
        time <- rep(x[3], length(variant))
        ans <- data.frame(variant, time)
    } else {
        ans <- NULL
    }
    return(ans)
}
plot_logihist <- function(place, dat)
{
    dat1a <- subset(dat, Country==place)
    dat1b <- do.call(rbind, apply(dat1a[, c("N_other", "N_focal", "time")], 1, convert_to_bernoulli))

    plot_colors <- c("data"="#888888", "freq"="black", "fit"="dodgerblue")
    logihist(dat1b$time, dat1b$variant, fillb=plot_colors["data"], 
                          colob=plot_colors["data"],
                          breaks = seq(min(dat1b$time), max(dat1b$time)+1) - 0.5,
                          ylab2="num observations") + ylab(paste0("freq(", variant, ")")) +
                  ggtitle(place) + theme_light() +
                  geom_point(data=dat1a, aes(x=time, y=Freq_focal), color=plot_colors["freq"], size=1)
}
plot_cases <- function(place, dat)
{
    dat1 <- subset(dat, Country==place)
    p <- ggplot(data = dat1) +
         geom_point(aes(x=time, y=Cases)) +
         labs(title=place) +
         theme_bw()
    p
}

if (variant == "omicron") {
    end_daynum <- 63
} else if (variant == "delta") {
    end_daynum <- 105
}
plots_variant_day <- lapply(places, function(x) {
                              plot_logihist(x, rename(dat3, time=daynum)) +
                              xlab("daynum") + geom_vline(xintercept=end_daynum) })
plots_case_day <- lapply(places, function(x) {
                              plot_cases(x, rename(dat3, time=daynum)) +
                              xlab("daynum") + geom_vline(xintercept=end_daynum) })

plots_day <- lapply(1:length(places),
                    function(i) plots_case_day[[i]] | plots_variant_day[[i]])
pdf(file=paste0("results/", variant, "-day.pdf", sep=""), width=8, height=3.4)
plots_day
dev.off()

#--------------------------------------------------
# Trim end of timeseries
#--------------------------------------------------

# Want most of the sigmoid shape, but not too far

# using end_daynum identified by plotting, above
dat4 <- dat3 %>%
        filter(daynum < end_daynum)

# avoid any countries without sequence data out to end_daynum
dropme <- dat4 %>%
          group_by(Country) %>%
          filter(N_total > 0) %>%
          summarize(m = max(daynum)) %>%
          filter(m < (end_daynum-2)) %>%
          pull(Country)
# (could be more careful about keeping countries only missing the last day or so,
#  but checking by eye, this seems fine)

dat4 <- dat4 %>%
        filter(!(Country %in% dropme))
places <- places[!(places %in% dropme)]

#--------------------------------------------------
# Aggregate by week
#--------------------------------------------------

# Plot the data as they will be handled by stan

dat5 <- dat4 %>%
        mutate(weeknum = daynum %/% 7) %>%
        group_by(weeknum, Country, first_day) %>%
        summarize(N_focal = sum(N_focal), N_other = sum(N_other),
                  N_total = sum(N_total), Cases = sum(Cases)) %>%
        ungroup() %>%
        mutate(Freq_focal = N_focal / N_total)


plots_variant_week <- lapply(places, function(x) {
                               plot_logihist(x, rename(dat5, time=weeknum)) +
                               xlab("weeknum") })
plots_case_week <- lapply(places, function(x) {
                               plot_cases(x, rename(dat5, time=weeknum)) +
                               xlab("weeknum") })

plots_week <- lapply(1:length(places),
                     function(i) plots_case_week[[i]] | plots_variant_week[[i]])
pdf(file=paste0("results/", variant, "-week.pdf", sep=""), width=8, height=3.4)
plots_week
dev.off()

#--------------------------------------------------
# Record data for use elsewhere
#--------------------------------------------------

write_csv(dat4, sprintf("data/%s-day-%s.csv", variant, filedate))   # used in fit_renewal2.R
write_csv(dat5, sprintf("data/%s-week-%s.csv", variant, filedate))  # could be used instead
