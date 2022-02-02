# schnellerGP

schnellerGP is an R package for fast, near-exact Gaussian Process sampling
in O(n log^2 n) time using the HODLRlib library of fast matrix operations.

# Installing in R

To install this package in R, begin by downloading the source copy of the package from the 
[releases folder](https://github.com/wesleyburr/schnellerGP/releases) and 
install it using type = "source". You will need several R packages, namely
Rcpp and RcppArmadillo, but should require no external dependencies. Note that
HODLRlib will perform much faster if you have access to OpenMP and have installed
the Intel MKL BLAS-replacement library (details in the **Performance Vignette** below).

The C++ code from HODLRlib and the underpinnings of the methods require a C++11
compiler. Note that if you are on a Mac, the clang compiler does not seem to reliably
compile the package - we are still debugging the Makefile to force that option.
No issues have been observed on Linux installations, and Windows binaries are available
at the releases folder on GitHub.

Once you're ready, you can install the package manually via the command:

    install.packages("schnellerGP-0.1-1.tar.gz", type = "source")
    
or by using the **Install Packages** GUI interface in RStudio (use the drop-down
to select *Package Archive File*, and then select the .tar.gz). Update the file 
name as necessary to whatever the copy you downloaded was set to for versioning.

# Getting Started

Refer to the [vignette on performance and installation](installation.html)
for details on installation and how to compile from source. 

More vignettes to come.

# Contributing?

Examine our [contribution guidelines](Contributing.md) for more.
