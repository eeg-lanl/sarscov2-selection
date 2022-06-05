require(tidyverse)
require(lubridate)
require(patchwork)
require(logihist)
require(coda)
require(latex2exp)

variant <- "omicron"
filedate <- "mar11H"
stanout <- mget(load(paste0("results/", variant, "-", filedate, "-ren-stanout.rda"),
                     envir=(tmp=new.env())), envir=tmp)

get_hpd <- function(x)
{
    ans <- HPDinterval(mcmc(x), prob=0.95)
    return(c(ans[1], median(ans), ans[2]))
}

#--------------------------------------------------
# Selective advantage in each country
#--------------------------------------------------
# and covariate weights, if any

get_hpd_cat <- function(x)
{
    ans0 <- median(x)
    ans1 <- HPDinterval(mcmc(x), prob=0.95)
    ans2 <- HPDinterval(mcmc(x), prob=0.50)
    ans <- data.frame(mid=ans0,
                      low1=ans1[1], high1=ans1[2],
                      low2=ans2[1], high2=ans2[2])
    return(ans)
}

# write a table of country s estimates
sel <- bind_rows(apply(stanout$chain[["s"]], 2, get_hpd_cat))
sel <- data.frame(Country = stanout$regions, sel)
write_csv(sel, file=file.path("results/", paste0(variant, "-", filedate, "-s.csv")))

# plot selective advantage estimates
plot_par_country <- function(parname, mname)
{
    X <- stanout$chain[[parname]]
    if (!is.na(mname))
    {
        mX <- stanout$chain[[mname]]
        X <- t(rbind(t(X), mX))
        cat_X <- bind_rows(apply(X, 2, get_hpd_cat)) %>%
                 add_column(label = c(stanout$regions, mname))
    } else {
        mname <- "mean"
        cat_X <- bind_rows(apply(X, 2, get_hpd_cat)) %>%
                 add_column(label = stanout$regions)

        cat_X <- cat_X %>%
                 add_row(!!! setNames(list(mean(cat_X$mid), NA, NA, NA, NA, mname), names(.)))
    }
    lev <- cat_X$label[order(cat_X$mid, decreasing=F)]
    lev <- rev(c(lev[lev!=mname], mname))

    cat_X <- cat_X %>% mutate(label = factor(label, levels=lev, ordered=T))

    X_plot <- ggplot(cat_X, aes(y=label, x=mid)) +
              geom_point(size=3) +
              geom_linerange(aes(xmin=low1, xmax=high1), size=0.5) +
              geom_linerange(aes(xmin=low2, xmax=high2), size=1.5) +
              theme_bw() +
              labs(title="", y="", x=parname)
    X_plot
}

s_cat_plot <- plot_par_country("s", "s_mean")
if (filedate == "mar11B")
    s_cat_plot <- plot_par_country("s", NA)

### Covariates ###

if (with(stanout$chain, exists("covariate_weights")))
{
    wts <- data.frame(wts = stanout$chain$covariate_weights[,1])
    wts_hpd <- get_hpd_cat(wts$wts)

    wts_plot <- ggplot(data = wts) +
                geom_density(aes(x = wts), size = 2) +
                geom_vline(xintercept = c(unlist(wts_hpd)), color="pink") +
                geom_vline(xintercept = 0, color="gray") +
                theme_bw()
}

pdf(file=paste0("results/", variant, "-", filedate, "-s.pdf", sep=""), width=5, height=8)
s_cat_plot
if (exists("wts_plot")) plot(wts_plot)
dev.off()

##############
# Timeseries #
##############
# cases, variant counts
# first assemble the summarized data, then plot in various ways

get_hpd_df <- function(x, prefix)
{
    ans <- data.frame(t(apply(x, 1, function(y) get_hpd(y))))
    colnames(ans) <- paste0(prefix, ".", c("low", "mid", "high"))
    return(ans)
}

#--------------------------------------------------
# Timeseries data: aggregated by week
#--------------------------------------------------

get_aggr_cases <- function(idx)
{
    nobs <- stanout$data_list$NumObs[stanout$regions[idx]]
    region <- stanout$regions[idx]
    scaling <- stanout$mean_cases$scaling[idx]
    cases.hat <- t(stanout$chain$HatCases[,1:nobs,idx]) * scaling
    cases.sim <- t(stanout$chain$SimCases[,1:nobs,idx]) * scaling

    ans <- data.frame(Country = region,
                      daynum = stanout$data_list$ObsTimes[1:nobs, region],
                      cases = stanout$data_list$Cases[1:nobs, region] * scaling,
                      get_hpd_df(cases.hat, "cases.hat"),
                      get_hpd_df(cases.sim, "cases.sim"))
    return(ans)
}
weekly_cases <- bind_rows(lapply(1:length(stanout$regions), get_aggr_cases))

get_aggr_variant <- function(idx)
{
    nobs <- stanout$data_list$NumObs[stanout$regions[idx]]
    region <- stanout$regions[idx]
    N_focal <- stanout$data_list$VariantSeq[1:nobs, region]
    N_total <- stanout$data_list$TotalSeq[1:nobs, region]
    freq.sim <- t(stanout$chain$SimVariantFreq[,1:nobs,idx])
    freq.hat <- t(stanout$chain$HatVariantFreq[,1:nobs,idx])

    ans <- data.frame(Country = region,
                      daynum = stanout$data_list$ObsTimes[1:nobs, region],
                      N_focal,
                      N_total,
                      freq = N_focal / N_total,
                      get_hpd_df(freq.sim, "freq.sim"),
                      get_hpd_df(freq.hat, "freq.hat"))
    return(ans)
}
weekly_variant <- bind_rows(lapply(1:length(stanout$regions), get_aggr_variant))

weekly_data <- tibble(full_join(weekly_cases, weekly_variant))
rm(weekly_cases, weekly_variant)

#--------------------------------------------------
# Timeseries data: by day
#--------------------------------------------------

get_day_cases <- function(idx)
{
    n <- stanout$data_list$N[stanout$regions[idx]]
    region <- stanout$regions[idx]
    scaling <- stanout$mean_cases$scaling[idx]
    cases.I <- t(stanout$chain$I[,,idx]) * scaling

    ans <- data.frame(Country = region,
                      daynum = 1:n,
                      cases = stanout$DayCases[[region]][1:n] * scaling,
                      get_hpd_df(cases.I, "cases.I"))
    return(ans)
}
daily_cases <- bind_rows(lapply(1:length(stanout$regions), get_day_cases))

get_day_variant <- function(idx)
{
    n <- stanout$data_list$N[stanout$regions[idx]]
    region <- stanout$regions[idx]
    N_total <- stanout$DayTotalSeq[[region]][1:n]
    N_focal <- stanout$DayVariantSeq[[region]][1:n]
    freq <- N_focal / N_total
    freq.pt <- t(stanout$chain$pt[,,idx])

    ans <- data.frame(Country = region,
                      daynum = 1:n,
                      N_focal,
                      N_total,
                      freq,
                      get_hpd_df(freq.pt, "freq.pt"))
    return(ans)
}
daily_variant <- bind_rows(lapply(1:length(stanout$regions), get_day_variant))

get_day_Rt <- function(idx)
{
    n <- stanout$data_list$N[stanout$regions[idx]] - 1
    region <- stanout$regions[idx]
    pt <- t(stanout$chain$pt[,1:n,idx])
    Rt_wt <- t(stanout$chain$Rt[,1:n,idx])
    fitness <- stanout$chain$fitness[,idx] 
    Rt <- ((1-pt) + pt*fitness)*Rt_wt

    ans <- data.frame(Country = region,
                      daynum = 1:n,
                      get_hpd_df(Rt_wt, "Rt_wt"),
                      get_hpd_df(Rt, "Rt"))
    return(ans)
}
daily_Rt <- bind_rows(lapply(1:length(stanout$regions), get_day_Rt))
  
daily_data <- tibble(daily_cases %>% full_join(daily_variant) %>% full_join(daily_Rt))
rm(daily_cases, daily_variant, daily_Rt)

#--------------------------------------------------
# Timeseries plots
#--------------------------------------------------

# (for logihist)
convert_to_bernoulli <- function(x)
{
    val <- c("old"=0, "new"=1)
    variant <- c(rep(val["old"], x[1]), rep(val["new"], x[2]))
    if (length(variant) > 0)
    {
        daynum <- rep(x[3], length(variant))
        ans <- data.frame(variant, daynum)
    } else {
        ans <- NULL
    }
    return(ans)
}

plot_country <- function(region, daily_data, weekly_data)
{
    dat_d <- daily_data %>% filter(Country == region)
    dat_w <- weekly_data %>% filter(Country == region)

    cls <- c("black", "dodgerblue")

    ### Cases ###

    # daily predictions and...
    p0 <- ggplot(data = dat_d, aes(x = daynum)) +
          geom_ribbon(aes(ymin = cases.I.low, ymax = cases.I.high), fill = cls[1], alpha = 0.3) +
          geom_line(aes(y = cases.I.mid), color = cls[1]) +
          labs(x = "day", y = "cases", title = region) +
          theme_bw()
          # scale_y_log10() + coord_cartesian(ylim=c(min(df$cases.hat.low) * 0.9, max(df$cases.hat.high) * 1.1))
    # ...daily data
    p1 <- p0 +
          geom_point(aes(y = cases), color = cls[1]) +
          labs(subtitle = "daily data, daily predictions")
    # ...weekly data [note the /7 hack]
    p2 <- p0 +
          geom_point(data = dat_w, aes(y = cases / 7), color = cls[1]) +
          labs(subtitle = "weekly data (/7), daily predictions")

    # weekly predictions and weekly data
    p3 <- ggplot(data = dat_w, aes(x = daynum)) +
          geom_ribbon(aes(ymin = cases.sim.low, ymax = cases.sim.high), fill = cls[1], alpha = 0.2) +
          geom_ribbon(aes(ymin = cases.hat.low, ymax = cases.hat.high), fill = cls[1], alpha = 0.3) +
          geom_line(aes(y = cases.hat.mid), color = cls[1]) +
          geom_point(aes(y = cases), color = cls[1]) +
          labs(x = "day", y = "cases", title = region) +
          labs(subtitle = "weekly data, weekly predictions") +
          theme_bw()

    ### Variant ###

    # daily data and daily predictions
    df2 <- do.call(rbind, apply(with(dat_d, data.frame(N_total-N_focal, N_focal, daynum)), 1, convert_to_bernoulli))
    p0 <- logihist(df2$daynum, df2$variant, fillb="#888888", colob="#888888",
                          breaks = seq(min(dat_d$daynum), max(dat_d$daynum)+1) - 0.5,
                          ylab2="num observations") +
          labs(x = "daynum", y = paste0("freq(", variant, ")"), title = region) +
          labs(subtitle = "daily data, daily predictions")
          # xlab("") + scale_x_continuous(labels = function(x) x+first_day) +
          # theme(axis.text.x=element_text(angle=25, hjust=1))
    p4 <- p0 +
          geom_ribbon(data=dat_d, aes(x=daynum, ymin=freq.pt.low, ymax=freq.pt.high), fill=cls[2], alpha=0.4, inherit.aes=F) +
          geom_line(data=dat_d, aes(x=daynum, y=freq.pt.mid), color=cls[2], inherit.aes=F) +
          geom_point(data=dat_d, aes(x=daynum, y=freq), color=cls[1], size=1.5, inherit.aes=F) +
          theme_bw()

    # weekly data and...
    df2 <- do.call(rbind, apply(with(dat_w, data.frame(N_total-N_focal, N_focal, daynum)), 1, convert_to_bernoulli))
    p0 <- logihist(df2$daynum, df2$variant, fillb="#888888", colob="#888888",
                          breaks = seq(min(dat_w$daynum), max(dat_w$daynum)+1) - 0.5,
                          ylab2="num observations") +
          labs(x = "daynum", y = paste0("freq(", variant, ")"), title = region)
    # ...weekly predictions
    p5 <- p0 +
          geom_ribbon(data=dat_w, aes(x=daynum, ymin=freq.sim.low, ymax=freq.sim.high), fill=cls[2], alpha=0.3, inherit.aes=F) +
          geom_ribbon(data=dat_w, aes(x=daynum, ymin=freq.hat.low, ymax=freq.hat.high), fill=cls[2], alpha=0.4, inherit.aes=F) +
          geom_line(data=dat_w, aes(x=daynum, y=freq.hat.mid), color=cls[2], inherit.aes=F) +
          geom_point(data=dat_w, aes(x=daynum, y=freq), color=cls[1], size=1.5, inherit.aes=F) +
          theme_bw() +
          labs(subtitle = "weekly data, weekly predictions")
    # ...daily predictions
    p6 <- p0 +
          geom_ribbon(data=dat_d, aes(x=daynum, ymin=freq.pt.low, ymax=freq.pt.high), fill=cls[2], alpha=0.4, inherit.aes=F) +
          geom_line(data=dat_d, aes(x=daynum, y=freq.pt.mid), color=cls[2], inherit.aes=F) +
          geom_point(data=dat_d, aes(x=daynum, y=freq), color=cls[1], size=1.5, inherit.aes=F) +
          theme_bw() +
          labs(subtitle = "weekly data (daily freqs), daily predictions")

    ### Rt ###

    # daily predictions
    p7 <- ggplot(data = dat_d, aes(x = daynum)) +
          geom_ribbon(aes(ymin=Rt_wt.low, ymax=Rt_wt.high), fill=cls[1], alpha=0.3) +
          geom_line(aes(y=Rt_wt.mid), color=cls[1]) +
          geom_ribbon(aes(ymin=Rt.low, ymax=Rt.high), fill=cls[2], alpha=0.3) +
          geom_line(aes(y=Rt.mid), color=cls[2]) +
          geom_hline(yintercept=1.0, color="gray10") +
          coord_cartesian(y=c(0, 3)) +
          labs(x = "daynum", y = TeX("$R_t$"), title = region) +
          labs(subtitle = "daily predictions") +
          theme_bw()

    # return(list(p1, p2, p3, p4, p6, p5, p7))
    return(list(p1, p3, p4, p5, p7))
}

all_plots <- lapply(stanout$regions, plot_country, daily_data, weekly_data)
all_plots_paged <- lapply(all_plots, wrap_plots, ncol=2)
pdf(file=paste0("results/", variant, "-", filedate, "-fit.pdf", sep=""), width=10, height=10)
all_plots_paged
dev.off()

quit()

#--------------------------------------------------
# Mean-squared error
#--------------------------------------------------

freq_mse <- function(country, daily_data, weekly_data)
{
    dat_d <- daily_data %>% filter(Country == country)
    dat_w <- weekly_data %>% filter(Country == country)

    # daily data and daily predictions
    day_mse <- mean((dat_d$freq - dat_d$freq.pt.mid)^2, na.rm=T)
    day_abs <- mean(abs(dat_d$freq - dat_d$freq.pt.mid), na.rm=T)

    # weekly data and weekly predictions
    week_mse <- mean((dat_w$freq - dat_w$freq.hat.mid)^2, na.rm=T)
    week_abs <- mean(abs(dat_w$freq - dat_w$freq.hat.mid), na.rm=T)

    return(tibble(country, day_mse, week_mse, day_abs, week_abs))
}

all_mse <- bind_rows(lapply(stanout$regions, freq_mse, daily_data, weekly_data))

colMeans(all_mse[,-1])
#     day_mse    week_mse     day_abs    week_abs
# 0.008882505 0.002235890 0.051891518 0.024588401

p1 <- ggplot(data = all_mse) + geom_histogram(aes(x=day_mse)) + theme_bw()
p2 <- ggplot(data = all_mse) + geom_histogram(aes(x=week_mse)) + theme_bw()
p1 + p2 + plot_annotation(title = "Mean-squared error of Omicron frequency")

layers <- list(theme_bw(),
               coord_cartesian(x = c(0, 0.125), y = c(0, 6)),
               labs(x = "avg abs error", y = "num countries"))
p1 <- ggplot(data = all_mse) +
      geom_histogram(aes(x=day_abs), binwidth = 0.004) +
      labs(subtitle = "daily data") +
      layers
p2 <- ggplot(data = all_mse) +
      geom_histogram(aes(x=week_abs), binwidth = 0.004) +
      labs(subtitle = "weekly data") +
      layers

pdf(file = "results/omicron_err.pdf", width=6, height=3)
p1 + p2 + plot_annotation(title = "average absolute error in Omicron frequency")
dev.off()

#--------------------------------------------------
# Example countries
#--------------------------------------------------

all_plots <- lapply(c("United Kingdom", "Switzerland", "India", "Brazil"), plot_country, daily_data, weekly_data)

plotme <- sapply(all_plots, function(x) x[c(2,4,5)])
plotme <- lapply(plotme, function(x) x + theme(plot.subtitle=element_blank()))
plotme[c(2,3,5,6,8,9,11,12)] <- lapply(plotme[c(2,3,5,6,8,9,11,12)], function(x) x + theme(plot.title=element_blank()))
plotme[4:12] <- lapply(plotme[4:12], function(x) x + theme(axis.title.y=element_blank()))
plotme[c(1:2,4:5,7:8,10:11)] <- lapply(plotme[c(1:2,4:5,7:8,10:11)], function(x) x + theme(axis.title.x=element_blank()))

plotme[c(1,4,7,10)] <- lapply(plotme[c(1,4,7,10)], function(x) x + scale_y_continuous(labels = scales::scientific))

wrap_plots(plotme, byrow=F)

pdf(file=paste0("results/", variant, "-", filedate, "-examples.pdf", sep=""), width=10, height=5)
wrap_plots(plotme, byrow=F)
dev.off()

#--------------------------------------------------
# NPI
#--------------------------------------------------
# not using this now

npi0 <- read_csv("data/government_response_index.csv")

npi1 <- npi0 %>%
        filter(country_name %in% sel.regions) %>%
        select(-c(1, country_code))

dat_npi <- npi1 %>%
           pivot_longer(cols = -c("country_name"), names_to="day1") %>%
           mutate(day = dmy(day1)) %>%
           filter(day >= first_day & day <= last_day)

plot_NPI <- function(region)
{
    start_day = start.dates[[region]]
    dat <- dat_npi %>%
           filter(country_name == region) %>%
           filter(day >= start_day)

    p <- ggplot(data = dat) +
         geom_line(aes(x=day, y=value), size=2) +
         xlab("day") + ylab("overall NPI") +
         coord_cartesian(y=c(0, 100)) +
         ggtitle(region) +
         theme_bw()
}

NPI_plot_list <- lapply(sel.regions, plot_NPI)
names(NPI_plot_list) <- sel.regions

stitch_plots <- function(region)
{
    p1 <- cases_plot_list[[region]] + ggtitle("")
    p2 <- freq_plot_list[[region]] + ggtitle("")
    p3 <- (NPI_plot_list[[region]] + ggtitle("")) / (Rt_plot_list[[region]] + ggtitle("")) + plot_layout(heights=c(1, 4))
    p <- (p1 | p2 | p3) + plot_annotation(title=region)
    p
}
