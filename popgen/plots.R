require(coda)
require(reshape2)
require(ggplot2)
require(gridExtra)
require(logihist)

#--------------------------------------------------
# Read in results (see fits.R)
#--------------------------------------------------

variant <- "R.1"
cases_so_far <- 2500

min_num <- 20

outname <- paste0(variant, "-", cases_so_far, "x", min_num)
load(paste0("stanfit-", outname, ".rda"))

fit_df <- as.data.frame(stanfit)

#--------------------------------------------------
# MCMC diagnostics
#--------------------------------------------------

plot(fit_df$lp__)
plot(fit_df$sgen_mean)

#--------------------------------------------------
# Make data frame for selection and migration
#--------------------------------------------------

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
            labs(title=paste(variant, cases_so_far, last_day, sep=", "), y="")
m_plot <- s_plot %+% cat_m

pdf(file=paste0("sm-", outname, ".pdf"), width=6, height=6)
s_plot + labs(x="s", subtitle="selection")
m_plot + labs(x="m", subtitle="migration")
dev.off()

#--------------------------------------------------
# Distribution of selection coefficients
#--------------------------------------------------

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

s_prior <- fit_df[,c("sgen_mean", "sgen_sd")]
s_prior <- na.omit(s_prior)
s_prior_samp <- apply(s_prior, 1, function(x) rnorm(n=10, mean=x[1], sd=x[2]))
s_prior_samp <- data.frame(s=as.vector(s_prior_samp))

p1 <- plot_ci(s_prior$sgen_mean)
p2 <- plot_ci(s_prior$sgen_sd)
p3 <- plot_ci(s_prior_samp$s)

more_layers <- list(geom_density(size=1.5),
                    theme_bw(),
                    ylab("probability density")
)

p1 <- p1 + more_layers + coord_cartesian(xlim=c(-1, 1.25)) + xlab("s mean")
p2 <- p2 + more_layers + coord_cartesian(xlim=c(0, 1)) + xlab("s std dev") + ylab("")
p3 <- p3 + more_layers + coord_cartesian(xlim=c(-1, 1.25)) + xlab("s") + ylab("")

pdf(file=paste0("selection-", outname, ".pdf", sep=""), width=8, height=2.25)
grid.arrange(p1, p2, p3, nrow=1, top=paste0(variant, ", ", cases_so_far, " cases"))
dev.off()

#--------------------------------------------------
# Data and predicted probabilities
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

# (for summarizing predicted freqs)
get_hpd <- function(fit_df)
{
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
                                   ylab2="num observations") + ylab(paste("freq(", variant, ")", sep="")) + 
                          xlab("") + scale_x_continuous(labels = function(x) x+first_day) +
                          ggtitle(place) + theme_light() +
                          theme(axis.text.x=element_text(angle=25, hjust=1)) +
                          geom_point(data=dat1a, aes(x=daynum, y=freq_new), color=plot_colors["freq"], size=0.9)
}

### Fit ###

p95 <- get_hpd(fit_df)
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

pdf(file=paste0("logihist-", outname, ".pdf", sep=""), width=4, height=3.4)
for (place in places)
{
    print(plot_data[[place]])
}
dev.off()
