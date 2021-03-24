library(reshape2)
library(countrycode)
library(tidyverse)
library(xlsx)

#GISAID metadata
d = read_lines("../data/metadata_2021-02-19_07-55.tsv.gz")
d = strsplit(d, "\t")
date = as.Date(sapply(d[2:length(d)], "[[", 5))
country = countrycode(sapply(d[2:length(d)], "[[", 7), origin = "country.name", destination = "iso3c")
#country = countrycode(country, origin = "iso3c", destination = "country.name")
host = sapply(d[2:length(d)], "[[", 15)
age = sapply(d[2:length(d)], "[[", 16)
sex = sapply(d[2:length(d)], "[[", 17)
pl = sapply(d[2:length(d)], "[[", 19)

dat = data.frame(date, country, host, age, sex, pl)
dat = dat[!is.na(dat$date),] 
upl = unique(dat$pl)

dat$target = dat$pl==upl[1]
ret = group_by(dat, country, date) %>%
  summarise(N=length(target), TEMP=sum(target))
names(ret)[4] = upl[1]

for(u in upl[2:length(upl)]){
  print(u)
  dat$target = dat$pl==u
  tmp = group_by(dat, country, date) %>%
    summarise(TEMP=sum(target))
  names(tmp)[3] = u
  ret = merge(ret, tmp)
}


#Hopkins covid epi data
dth = read.csv("./COVID-19-master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv", header=T)
dth1 = aggregate(dth[,5:410], by=list(dth$Country.Region), FUN=sum)
rownames(dth1) = dth1$Group.1
dth1 = as.data.frame(t(apply(dth1, 1, function(x)(diff(c(0, as.numeric(x[2:407])))))))
names(dth1) = names(dth)[5:410]
dth = melt(t(dth1))
dth$date = as.Date(dth$Var1, format="X%m.%d.%y")
dth$Var1 = NULL
dth$Var2 = countrycode(as.character(dth$Var2), origin = "country.name", destination = "iso3c")
names(dth)[1:2] = c("country", "deaths")

cas = read.csv("./COVID-19-master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv", header=T)
cas1 = aggregate(cas[,5:410], by=list(cas$Country.Region), FUN=sum)
rownames(cas1) = cas1$Group.1
cas1 = as.data.frame(t(apply(cas1, 1, function(x)(diff(c(0, as.numeric(x[2:407])))))))
names(cas1) = names(cas)[5:410]
cas = melt(t(cas1))
cas$date = as.Date(cas$Var1, format="X%m.%d.%y")
cas$Var1 = NULL
cas$Var2 = countrycode(as.character(cas$Var2), origin = "country.name", destination = "iso3c")
names(cas)[1:2] = c("country", "cases")

rec = read.csv("./COVID-19-master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv", header=T)
rec1 = aggregate(rec[,5:410], by=list(rec$Country.Region), FUN=sum)
rownames(rec1) = rec1$Group.1
rec1 = as.data.frame(t(apply(rec1, 1, function(x)(diff(c(0, as.numeric(x[2:407])))))))
names(rec1) = names(rec)[5:410]
rec = melt(t(rec1))
rec$date = as.Date(rec$Var1, format="X%m.%d.%y")
rec$Var1 = NULL
rec$Var2 = countrycode(as.character(rec$Var2), origin = "country.name", destination = "iso3c")
names(rec)[1:2] = c("country", "recoveries")

dat2 = merge(merge(dth, cas, all=T), rec, all=T)

ret$date = as.Date(ret$date)
dat2$date = as.Date(dat2$date)
#dat$country = countrycode(dat$country, origin = "iso3c", destination = "country.name")
#dat2$country = countrycode(dat2$country, origin = "iso3c", destination = "country.name")

d = merge(ret, dat2, all=T)
d$N[is.na(d$N)] = 0
d$N_mt[is.na(d$N_mt)] = 0
d$country = countrycode(d$country, origin = "iso3c", destination = "country.name")
write.table(d, file="seq-pl-all_epi-26022021", sep=",", row.names = F)

#===========================================
library(tidyverse)
library(reshape2)

d = read.csv("seq-pl-all_epi-26022021")
d = pivot_longer(d, 4:886)
d$pr = d$value/d$N
d = d[d$pr>0,]
agg = group_by(d, country) %>% summarise(nms=name[])

