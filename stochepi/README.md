# Simulate a 2-variant SARS-CoV-2 model and fit to data with SMC

## Compiling the code

The code is written in C++ version C++17. It is tested with the GCC compiler
version 9.3.0. Compile the code using `make`.

## Running the program

Run the program `estavoir` from the command line. The command line option `-h`
shows a list of options.

## Parameter and data files

Parameter specifications have to be put in a file (e.g. `parameters.txt`), and
the filename has to be specified with the `-p` option on the command line.
The parameter file contains the initial parameter guess, lower and upper bounds,
and the standard deviation of the random walk. These values have to be separated
by spaces. Empty lines and lines starting with `#` are ignored. Example:

```
## specify parameters as follows:
## name initial_guess sigma_rw lower_bound upper_bound

beta 0.1 1e-3 0 inf
gamma 1.2 0 -1 1
```
the keywords `inf` and `-inf` are used to indicate unbounded parameters.
When the standard deviation for the random walk is `0`, the parameter is fixed.

Data has to be given in a tab-separated-values file (e.g.`data.tsv`), and
the file name has to be given with the `-d` option.
The data file must be organized as follows:
```
ID  time  event  deaths  deaths_cc var_seq var_seq_cc  total_seq total_seq_cc
```
where `ID` is the name of the region, `time` is the time of the observation,
`deaths` is the number of deaths between this and the previous observation,
`deaths_cc` is a code determining the censoring of deaths, `var_seq` is the
number of variant sequences, `total_seq` is the total number of sampled sequences
and `event` is a string that signifies an event at the observation time.

The censoring codes are `0` for uncensored, `1` for left-censored, `2` for
right-censored and `3` for missing data. The `event` string has to be
`[RESET_CASES]`, which makes sure that the accumulator variables are set to
zero after an observation.

## Interpreting results with jupyter notebooks

In the `notebooks` folder, the following jupyter notebook can be found:

* `ParseData.ipynb` prepare input data files for SMC.
* `FigSMCFitSingle.ipynb` check the output of SMC for a single variant and region, and look at diagnostics and parameter estimates.
* `ExtractProfileLik.ipynb` extract the log-likelihood from a batch of SMC fits to construct a profile likelihood.
* `FigProfLik.ipynb` Combine profile likelihoods from mutliple regions and variants into a single figure.
* `ExtractProfLikeTmax.ipynb` again, extract log-likelihoods, but now for multiple time horizons.
* `FigSMCFitTmaxAnalysisMultiRegion.ipynb` make a figure with the SMC fit, combined with the profile likelihoods for multiple time horizons.
* `ReffAnalysis.ipynb` Compute the effective reproduction number as a function of time. This notebook can also be used to compute the external force of infection for the migration model
* `CompareModels.ipynb` Analysis of how different models (population genetics vs. epidemic model) lead to different estimates of the selective advantage $s$.
