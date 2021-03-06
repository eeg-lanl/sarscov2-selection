#--------------------------------------------------
# Read in results (see fits.R)
#--------------------------------------------------

### B.1.1.7

load("fits-b117.rda")
label <- "B.1.1.7"

### D614G

load("fits-614.rda")
label <- "D614G"

#--------------------------------------------------
# Make data frame for selection and migration
#--------------------------------------------------

require(reshape)

fit_df <- as.data.frame(fitN)

i <- grep("sgen\\[", names(fit_df))
sgen <- fit_df[,i]
names(sgen) <- places
sgen <- melt(sgen)
names(sgen) <- c("Country", "s")

i <- grep("mgen\\[", names(fit_df))
mgen <- fit_df[,i]
names(mgen) <- places
mgen <- melt(mgen)
names(mgen) <- c("Country", "m")

all(sgen$Country == mgen$Country)
fit_sm <- cbind(sgen, m=mgen$m)

#--------------------------------------------------
# Caterpillars for selection and migration
#--------------------------------------------------

require(coda)
require(ggplot2)

get_hpd <- function(x)
{
    ans0 <- median(x)
    ans1 <- HPDinterval(mcmc(x), prob=0.95)
    ans2 <- HPDinterval(mcmc(x), prob=0.90)
    ans <- c(ans0, ans1, ans2)
    names(ans) <- c("mid", "low1", "high1", "low2", "high2")
    return(ans)
}

cat_s <- data.frame(do.call("rbind", by(fit_sm$s, fit_sm$Country, get_hpd)))
cat_s$country <- rownames(cat_s)
rownames(cat_s) <- NULL

cat_m <- data.frame(do.call("rbind", by(fit_sm$m, fit_sm$Country, get_hpd)))
cat_m$country <- rownames(cat_m)
rownames(cat_m) <- NULL

cat_m$country <- cat_s$country <- factor(cat_s$country, levels=cat_s$country[order(cat_s$mid, decreasing=F)], ordered=T)

s_plot <- ggplot(cat_s, aes(y=country, x=mid)) + geom_point(size=3) + 
            geom_linerange(aes(xmin=low1, xmax=high1), size=0.5) + 
            geom_linerange(aes(xmin=low2, xmax=high2), size=1.5) + 
            labs(title=label, y="")
m_plot <- s_plot %+% cat_m

pdf(file=paste("caterpillar-", label, ".pdf", sep=""), width=6, height=6)
s_plot + labs(x="s", subtitle="selection")
m_plot + labs(x="m", subtitle="migration")
dev.off()

#--------------------------------------------------
# Correlation between selection and migration
#--------------------------------------------------

require(ggplot2)

pl <- ggplot(data=fit_sm, aes(x=s, y=m)) + geom_bin2d(bins=50, show.legend=F) +
        scale_fill_continuous(type="viridis") + facet_wrap(~Country) +
        labs(title=label, subtitle="parameter correlations", x="selection, s", y="migration, m")

pdf(file=paste("correlation-", label, ".pdf", sep=""), width=12, height=10)
if (label == "B.1.1.7")
{
    pl + ylim(-0.01, 0.08) + xlim(-1, 1)
} else {
    pl
}
dev.off()

#--------------------------------------------------
# Distribution of selection coefficients
#--------------------------------------------------

require(ggplot2)
require(cowplot)
require(coda)

# medium.com/@jireh/a-clever-use-of-ggplot-internals-bbb168133909
plot_ci <- function(x)
{
    gg_density <- ggplot(data=data.frame(x), aes(x=x)) + geom_density()
    ci <- HPDinterval(mcmc(x), prob=0.9)
    mid <- median(x)

    build_object <- ggplot_build(gg_density)
    x_dens <- build_object$data[[1]]$x
    y_dens <- build_object$data[[1]]$y
    index_left <- min(which(x_dens >= ci[1]))
    index_right <- max(which(x_dens <= ci[2]))
    index_mid <- which.min(abs(x_dens - mid))
    y_mid <- y_dens[index_mid]

    gg_density + 
        geom_area(
            data=data.frame(
              x=x_dens[index_left:index_right],
              y=y_dens[index_left:index_right]), 
            aes(x=x,y=y), fill="#CCCCCC") +
        geom_segment(
            aes(x=mid, xend=mid, y=0, yend=y_mid), 
            color="#555555")
}

### The results for the top-level distribution on selection

s_prior <- fit_df[,c("sgen_prior_mean", "sgen_prior_sd")]
s_prior_samp <- apply(s_prior, 1, function(x) rnorm(n=10, mean=x[1], sd=x[2]))
s_prior_samp <- data.frame(s=as.vector(s_prior_samp))

p1 <- plot_ci(s_prior$sgen_prior_mean)
p2 <- plot_ci(s_prior$sgen_prior_sd)
p3 <- plot_ci(s_prior_samp$s)

more_layers <- list(geom_density(size=1.5),
                    theme_bw(),
                    ylab("probability density")
)

p1 <- p1 + more_layers + xlim(0, 0.5) + xlab("s mean")
p2 <- p2 + more_layers + xlim(0, 0.5) + xlab("s std dev") + ylab("")
p3 <- p3 + more_layers + xlim(-0.75, 1.25) + xlab("s") + ylab("")

lab <- list("D614G"=c("A", "B", "C"), "B.1.1.7"=c("D", "E", "F"))

pdf(file=paste("selection-", label, ".pdf", sep=""), width=8, height=2)
plot_grid(p1, p2, p3, nrow=1, labels=lab[[label]])
dev.off()

### For the focal countries

require(bayesplot)
color_scheme_set("gray")

p <- c("United Kingdom", "Netherlands")
i <- match(p, places)
s_places <- fit_df[, paste("sgen[", i, "]", sep="")]
names(s_places) <- p

pdf(file=paste("selection-", label, "-2.pdf", sep=""), width=4.5, height=4)
mcmc_areas(s_places, prob=0.9) + 
    labs(title=label, subtitle="selection coefficient") +
    yaxis_ticks(on=F) + yaxis_text(on=T) + xlim(c(0, 1))
dev.off()

#--------------------------------------------------
# Data and predicted probabilities
#--------------------------------------------------

require(logihist)
require(coda)

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

# (for summarizing predicted freqs)
get_hpd <- function(fit_stan)
{
    fit_df <- as.data.frame(fit_stan)
    i_theta <- grep("theta", names(fit_df))
    ans <- data.frame(t(apply(fit_df[,i_theta], 2, function(x) c(median(x), HPDinterval(mcmc(x))))))
    names(ans) <- c("med", "low", "high")
    return(ans)
}

### Data ###

plot_colors <- c("data"="#888888", "freq"="black", "fit"="dodgerblue")

plot_data <- list()
for (place in places)
{
    dat1a <- subset(dat, country==place)
    dat1a <- subset(dat1a, num_old > 0 | num_new > 0)
    dat1b <- do.call(rbind, apply(dat1a[, c("num_old", "num_new", "daynum")], 1, convert_to_bernoulli))

    plot_data[[place]] <- logihist(dat1b$daynum, dat1b$variant, fillb=plot_colors["data"], 
                                   colob=plot_colors["data"], intervalo=1, 
                                   ylab2="num observations") + ylab(paste("freq(", label, ")", sep="")) + 
                          xlab("") + scale_x_continuous(labels = function(x) x+first_day) +
                          ggtitle(place) + theme_light() +
                          theme(axis.text.x=element_text(angle=25, hjust=1)) +
                          geom_point(data=dat1a, aes(x=daynum, y=freq_new), color=plot_colors["freq"], size=0.9)
}

### Fit

p95 <- get_hpd(fitN)
predicted <- data.frame(daynum=dat$daynum, p95, country=places[dat$id])

for (place in places)
{
    plot_data[[place]] <- plot_data[[place]] + 
                            geom_ribbon(data=subset(predicted, country==place), 
                                        aes(x=daynum, ymin=low, ymax=high), 
                                        fill=plot_colors["fit"], alpha=0.75, inherit.aes=F) +
                            geom_line(data=subset(predicted, country==place), 
                                      aes(x=daynum, y=med), color=plot_colors["fit"])
}

pdf(file=paste("logihist-", label, ".pdf", sep=""), width=4, height=3.4)
for (place in places)
{
    print(plot_data[[place]])
}
dev.off()

#--------------------------------------------------
# Compare s for both mutations in each country
#--------------------------------------------------

require(ggplot2)

# run code above...
s_b117 <- cat_s
s_d614g <- cat_s
dat <- merge(s_d614g, s_b117, by="country")

p <- ggplot(data=dat) + theme_bw() +
     geom_abline(slope=1, intercept=0, color="gray") +
     geom_pointrange(aes(x=mid.x, y=mid.y, xmin=low1.x, xmax=high1.x)) +
     geom_pointrange(aes(x=mid.x, y=mid.y, ymin=low1.y, ymax=high1.y)) +
     xlab("s for D614G") + ylab("s for B.1.1.7") +
     coord_cartesian(xlim=c(-0.25, 0.75), ylim=c(-0.25, 0.75))

pdf("s-vs-s.pdf", width=5, height=5)
p
dev.off()
