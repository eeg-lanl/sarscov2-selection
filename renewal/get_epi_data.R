## extract COVID epi data from github 
## working folder is assumed to be equal to the project root

library(countrycode)
library(reshape2)
library(zoo)
library(tidyverse)

download.file("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv",
                 destfile = "data/cases.csv")
download.file("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv",
                 destfile = "data/deaths.csv")
download.file("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv",
                 destfile = "data/recovery.csv")
cas = read.csv("data/cases.csv")
dth = read.csv("data/deaths.csv")
rec = read.csv("data/recovery.csv")
n.days = length(names(cas)) - 4

#Hopkins covid epi data
#dth = read.csv("../../COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv", header=T)
dth1 = aggregate(dth[,5:n.days], by=list(dth$Country.Region), FUN=sum) #these ranges need to be adjusted depending on temporal range of data
rownames(dth1) = dth1$Group.1
dth1 = as.data.frame(t(apply(dth1, 1, function(x)(diff(c(0, as.numeric(x[2:(n.days-3)]))))))) #these ranges need to be adjusted depending on temporal range of data
names(dth1) = names(dth)[5:n.days]#these ranges need to be adjusted depending on temporal range of data
dth = melt(t(dth1))
dth$date = as.Date(dth$Var1, format="X%m.%d.%y")
dth$Var1 = NULL
dth$Var2 = countrycode(as.character(dth$Var2), origin = "country.name", destination = "iso3c")
names(dth)[1:2] = c("country", "deaths")

#cas = read.csv("../../COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv", header=T)
cas1 = aggregate(cas[,5:n.days], by=list(cas$Country.Region), FUN=sum)#these ranges need to be adjusted depending on temporal range of data
rownames(cas1) = cas1$Group.1
cas1 = as.data.frame(t(apply(cas1, 1, function(x)(diff(c(0, as.numeric(x[2:(n.days-3)])))))))#these ranges need to be adjusted depending on temporal range of data
names(cas1) = names(cas)[5:n.days]#these ranges need to be adjusted depending on temporal range of data
cas = melt(t(cas1))
cas$date = as.Date(cas$Var1, format="X%m.%d.%y")
cas$Var1 = NULL
cas$Var2 = countrycode(as.character(cas$Var2), origin = "country.name", destination = "iso3c")
names(cas)[1:2] = c("country", "cases")

#rec = read.csv("../../COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv", header=T)
rec1 = aggregate(rec[,5:n.days], by=list(rec$Country.Region), FUN=sum)#these ranges need to be adjusted depending on temporal range of data
rownames(rec1) = rec1$Group.1
rec1 = as.data.frame(t(apply(rec1, 1, function(x)(diff(c(0, as.numeric(x[2:(n.days-3)])))))))#these ranges need to be adjusted depending on temporal range of data
names(rec1) = names(rec)[5:n.days]#these ranges need to be adjusted depending on temporal range of data
rec = melt(t(rec1))
rec$date = as.Date(rec$Var1, format="X%m.%d.%y")
rec$Var1 = NULL
rec$Var2 = countrycode(as.character(rec$Var2), origin = "country.name", destination = "iso3c")
names(rec)[1:2] = c("country", "recoveries")

dat2 = merge(merge(dth, cas, all=T), rec, all=T)
#ret = NULL
#ret$date = as.Date(ret$date)
dat2$date = as.Date(dat2$date)
#dat$country = countrycode(dat$country, origin = "iso3c", destination = "country.name")
#dat2$country = countrycode(dat2$country, origin = "iso3c", destination = "country.name")
#names(ret)[1:2] = c("region", "country")
#d = merge(ret, dat2, all=T)
#d$N[is.na(d$N)] = 0
#d$N_mt[is.na(d$N_mt)] = 0
dat2$country = countrycode(dat2$country, origin = "iso3c", destination = "country.name")
dat2 = dat2[!duplicated(dat2),]
dat2 = dat2[-which(is.na(dat2$country)),]

#rm7_cases = rm3_cases = rm7_deaths = rm3_deaths = NULL
dat3 = NULL
for(u in unique(dat2$country)){
  print(u)
  tmp = dat2[dat2$country==u,]
  tmp = tmp[!duplicated(tmp),]
  tmp = tmp[order(tmp$date, decreasing = T),]
  tmp$rm7_cases = rollapply(tmp$cases, 7, mean, partial=T)
  tmp$rm3_cases = rollapply(tmp$cases, 3, mean, partial=T)
  n <- length(tmp$cases)
  ks <- ksmooth(1:n, tmp$rm7_cases, kernel="normal", bandwidth = 3)
  tmp$ks_cases = ks$y
  tmp$rm7_deaths = rollapply(tmp$deaths, 7, mean, partial=T)
  tmp$rm3_deaths = rollapply(tmp$deaths, 3, mean, partial=T)
  dat3 = rbind(dat3, tmp)
}
write.table(dat3, file="data/epi-data.csv", sep=",", row.names = F)
