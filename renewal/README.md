# Renewal model

This code is for a deterministic model that uses variant counts and case counts to estimate the selective advantage of an emerging variant.

It requires [R](https://www.r-project.org/) and [Stan](https://mc-stan.org/).

Prepare the data:
  * `get_epi_data.R`
  * `parse_omicron_delta.R`
  * `data_renewal.R`

Fit the model:
  * `focalvariant-renewal-aggr.stan`: 
  * `fit_renewal.R`: 

Visualize the results:
  * `plot_renewal.R`: 
