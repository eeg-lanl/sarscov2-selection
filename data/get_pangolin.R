library(tidyverse)
library(lubridate)
library(countrycode)

# from gisaid
dat <- read_csv("metadata_2021_05_20.csv")
problems(dat)
# (warnings are because the Variant column didn't parse well -- that's okay)
dat <- select(dat, Name=`Virus name`, Date=`Collection date`, Location=Location, Pango=`Pango lineage`)

# check pango assignments
filter(dat, is.na(Pango)) # only a few
filter(dat, Pango == "None") # a lot
dat <- replace_na(dat, list(Pango="None"))

# country code
loc <- str_split(dat$Location, "/")
country <- sapply(loc, function(x) x[2])
country <- str_trim(country)
cc <- countrycode(country, origin="country.name", destination="iso3c")
cc <- countrycode(cc, destination="country.name", origin="iso3c")
dat <- add_column(dat, Country=cc)

# date
dat <- mutate(dat, Day=ymd(Date))

# keep just what we need
dat1 <- dat %>% drop_na(Country, Day) %>% filter(Day > ymd("2019-11-30"))
dat1 <- select(dat1, Country, Day, Pango)

# combine all counts for each country on each day
dat2 <- pivot_wider(dat1, names_from=Pango, values_from=Pango, values_fn=list(Pango=length), values_fill=0)
dat2 <- mutate(dat2, Total=rowSums(across(where(is.numeric))))

# alphabetize the variants
dat2 <- dat2 %>% select(sort(colnames(.))) %>% select(Country, Day, Total, everything())

# sort rows by Day, then by Country
dat2 <- arrange(dat2, Day, Country)

# record
write_csv(dat2, file="pango_2021-05-20.csv")
