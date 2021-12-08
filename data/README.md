# Data

We're not allowed to redistribute the actual data, so here are the steps:

  * Download the genetic variant data in the form of the nextstrain metadata file from [GISAID](https://www.gisaid.org)
  * Get the case count data compiled by the [JHU CSSE team](https://github.com/CSSEGISandData/COVID-19)
  * Use `merge_pangolin_all.R` to arrange the data

The resulting file will have one row per day per country, with columns for the numbers of cases by genetic variant, total number of sequences, deaths, recoveries, and cases.

For the popgen model, which uses only sequence counts, `get_pangolin.R` will arrange the data.
The resulting file will have one row per day per country, with one column for each pangolin variant.
