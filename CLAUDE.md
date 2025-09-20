# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Architecture

**libpulsatile** is a C++ library with R package bindings for Bayesian pulsatile hormone modeling. The project has a dual-build system:

1. **Standalone C++ library** (root directory) - For development and testing C++ components
2. **R package** (`R-package/` directory) - Production interface using Rcpp/RcppArmadillo

### Core Components
- **Data Structures** (`include/bp_datastructures/`): Patient, Population, Chains, PulseEstimates
- **MCMC Algorithms** (`include/bp_mcmc/`): Metropolis-Hastings, Birth-Death process, Proposal variance adaptation
- **Single Subject Model** (`include/bpmod_singlesubject/`): Main modeling interface
- **R Interface** (`R-package/src/singlesubject.cpp`): Rcpp exports for R integration

### Symbolic Links
The project uses symbolic links between the main library and R package. Windows users must clone with `git clone -c core.symlinks=true` and have appropriate permissions.

## Dependencies
- **R packages**: Rcpp, RcppArmadillo, devtools, testthat, RInside
- **System libraries**: gfortran (required for RcppArmadillo)
- Install R dependencies: `Rscript -e "install.packages(c('Rcpp', 'RcppArmadillo', 'devtools', 'testthat', 'RInside'))"`

## Build and Test Commands
- **Build library**: `make`
- **Build R package**: `cd R-package && make`
- **Run C++ tests**: `bin/tests` or single test: `bin/tests [test-name]`
- **Run R tests**: `cd R-package && Rscript -e 'devtools::test()'`
- **Check R package**: `cd R-package && Rscript -e 'devtools::check()'`
- **Install R package**: `cd R-package && Rscript -e 'devtools::install()'`
- **Generate Rcpp exports**: `cd R-package && Rscript -e "Rcpp::compileAttributes()"`

## Code Style Guidelines
- **Headers**: Use include guards: `GUARD_filename_h`
- **Naming**: ClassNames, function_names, local_variables, CONSTANTS
- **Formatting**: 2-space indentation, 80-character line length
- **C++ Standards**: C++11, use Rcpp/RcppArmadillo for R integration
- **R Style**: Follow tidyverse style guide, use `<-` for assignment
- **Imports**: Group imports by library (Rcpp, std, project headers)
- **Documentation**: Roxygen2 for R functions, header comments for C++
- **Error Handling**: `Rcpp::stop()` for C++ errors, `stopifnot()` for R
- **Tests**: Catch framework for C++, testthat for R

## Development Workflow
- The main algorithm is in `singlesubject_()` Rcpp export function
- C++ development should follow the existing modular structure in `include/` directories
- When adding new functionality, add corresponding tests in both `tests/` (C++) and `R-package/tests/` (R)
- Use `make clean && make` to ensure clean builds after header changes

## CI/CD Pipeline
- **GitHub Actions** workflows in `.github/workflows/`:
  - `ci.yml`: Builds and tests on Ubuntu and macOS with R 4.3.0
  - `coverage.yml`: Generates coverage reports using `covr` package
- Replaces legacy Travis CI configuration

## Common Issues
- If compile fails with "library 'gfortran' not found", install gfortran for your platform
- For macOS with Homebrew: `brew install gfortran` then ensure R can find it: `sudo ln -sf /opt/homebrew/bin/gfortran /usr/local/bin/gfortran`
- For macOS general guidance, see: https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/
- Use `-DNORINSIDE` flag if you want to skip RInside dependency
- Windows users: Ensure symbolic links are properly created during clone
- ARM64 Mac compatibility: Fixed Catch framework assembly instructions for Apple Silicon