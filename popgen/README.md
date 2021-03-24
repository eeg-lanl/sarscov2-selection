# Population genetics model

This code is for the deterministic population genetics model with selection and migration.
It requires [R](https://www.r-project.org/) and [Stan](https://mc-stan.org/).

  * `fits.R`: Choose which variant and then run this to fit the model across all countries.
  * `selection_migration.stan`: The actual model, called by `fits.R`.
  * `plots.R`: Visualize the results.
