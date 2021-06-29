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
The parameter fie contains the initial parameter guess, lower and upper bounds,
and the standard deviation of the random walk. These values have to be seperated
by spaces. Empty lines and lines starting with `#` are ignored. Example:

```
## specify parameters as folows:
## name initial_guess sigma_rw lower_bound upper_bound

beta 0.1 1e-3 0 inf
gamma 1.2 0 -1 1
```
the keywords `inf` and `-inf` are used to indicate unbounded parameters.
When the standard deviation for the random walk is `0`, the parameter is fixed.

*TODO: data file*


## Interpreting results with jupyter notebooks

*TODO*
