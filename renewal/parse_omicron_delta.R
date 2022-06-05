### input:  gisaid metadata file
### output: file with counts per country per day,
###         for "Omicron" (B.1.1.529) or "Delta" (B.1.617.2+AY.x) vs "Other"
###         can cut back to an earlier sampling date

library(archive)
library(tidyverse)
library(lubridate)
library(patchwork)

# full metadata
dat0 <- read_tsv(archive_read(archive="data/metadata_tsv_2022_03_11.tar.xz", file="metadata.tsv"),
                 col_types = cols(
                    `Virus name` = col_character(),
                    Type = col_character(),
                    `Accession ID` = col_character(),
                    `Collection date` = col_character(),
                    Location = col_character(),
                    `Additional location information` = col_character(),
                    `Sequence length` = col_double(),
                    Host = col_character(),
                    `Patient age` = col_character(),
                    Gender = col_character(),
                    Clade = col_character(),
                    `Pango lineage` = col_character(),
                    `Pangolin version` = col_date(format = ""),
                    Variant = col_character(),
                    `AA Substitutions` = col_character(),
                    `Submission date` = col_date(format = ""),
                    `Is reference?` = col_logical(),
                    `Is complete?` = col_logical(),
                    `Is high coverage?` = col_logical(),
                    `Is low coverage?` = col_logical(),
                    `N-Content` = col_double(),
                    `GC-Content` = col_double()
                    ))

# drop rows that nextstrain says are bad
ex <- read_csv("https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/exclude.txt", comment="#", col_names=F)
ex <- ex %>% mutate(`Virus name` = paste0("hCoV-19/", X1))
dat1 <- anti_join(dat0, ex)

# just the columns we need
dat1 <- dat1 %>%
        select(Day = `Collection date`, Location, Pango = `Pango lineage`, VariantString = Variant) %>%
        mutate(Day = ymd(Day))

# drop rows with useless dates
dat1 <- dat1 %>% filter(!is.na(Day)) %>%
                 filter(Day > ymd("2019-11-30"))

# drop rows that lack both Pango and VariantString
dat1 <- dat1 %>% filter(!((is.na(Pango) | Pango == "None") & is.na(VariantString)))

# convert to omicron vs delta vs other
# (if VariantString and Pango disagree, trust Pango)
dat1 <- dat1 %>%
        mutate(Variant = case_when(
                         Pango == "B.1.1.529" ~ "Omicron",
                         str_starts(Pango, fixed("BA.")) ~ "Omicron",
                         str_detect(VariantString, "Omicron") ~ "Omicron",
                         Pango == "B.1.617.2" ~ "Delta",
                         str_starts(Pango, fixed("AY.")) ~ "Delta",
                         str_detect(VariantString, "Delta") ~ "Delta",
                         TRUE ~ "Other"
                         ))

# get country
dat1 <- dat1 %>%
        mutate(Country = str_split_fixed(Location, fixed(" / "), 3)[,2])

# for omicron
dat3 <- dat1 %>%
        mutate(Variant = replace(Variant, Variant != "Omicron", "Other"))
# first appearance of omicron
dat3 %>% filter(Variant == "Omicron") %>% arrange(Day)
datO1 <- dat3 %>% filter(Day >= ymd("2021-11-01"))

# for delta
dat3 <- dat1 %>%
        mutate(Variant = replace(Variant, Variant != "Delta", "Other"))
# first appearance of delta
dat3 %>% filter(Variant == "Delta") %>% arrange(Day) %>% print(n=30)
# first 2020-03-12, gap until 2020-06-26, first in India 2020-10-23
# but first sequence wasn't submitted until 2021-03-18
datD1 <- dat3 %>% filter(Day >= ymd("2021-01-01"))

# per country per variant per day
lump_data <- function(dat)
{
    dat <- dat %>%
           group_by(Day, Country, Variant) %>%
           summarize(Count = n()) %>%
           ungroup()
    dat <- dat %>%
           pivot_wider(names_from = Variant, values_from = Count, values_fill = 0) %>%
           select(any_of(c("Day", "Country", "Delta", "Omicron", "Other"))) %>%
           arrange(Day, Country)
    return(dat)
}

datO2 <- lump_data(datO1)
datD2 <- lump_data(datD1)

# write the parsed data file
write_csv(datO2, "data/omicron-2022-03-11.csv")
write_csv(datD2, "data/delta-2022-03-11.csv")

# visual check
p0 <- ggplot(data=datD2, aes(x=Day)) + theme_bw()
p1 <- p0 + geom_col(aes(y=Delta))
p2 <- p0 + geom_col(aes(y=Other))
p_delta <- p1 + p2 + plot_layout(ncol=1)
p0 <- ggplot(data=datO2, aes(x=Day)) + theme_bw()
p1 <- p0 + geom_col(aes(y=Omicron))
p2 <- p0 + geom_col(aes(y=Other))
p_omicron <- p1 + p2 + plot_layout(ncol=1)
