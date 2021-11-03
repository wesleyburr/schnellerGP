# schnellerGP

schnellerGP is an R package for fast, near-exact Gaussian Process sampling
in O(n log^2 n) time using H and H^2 matrix decompositions. 

# Installing in R

To install this package in R, begin by downloading the source copy of the package from the 
[releases folder](https://github.com/wesleyburr/schnellerGP/releases) and 
install it using type = "source". To install from source, you will require several
libraries for your system: LAPACK, BLAS, gFortran, NetCDF and Cairo. Some attempt
has been made at describing how to install these for common architectures
in the **Installation Vignette**, below. 

The package has no external dependencies, but requires several libraries
that are built-in via source (C++ libraries for H2 (h2lib) and HODLR (hodlrlib)).
This requires (for a source compilation) a C++ compiler and a number of external
libraries, as mentioned above. 

Once you're setup, you can install the package manually via the command:

    install.packages("schnellerGP-0.1-0.tar.gz", type = "source")
    
or by using the **Install Packages** GUI interface in RStudio (use the drop-down
to select *Package Archive File*, and then select the .tar.gz). Update the file 
name as necessary to whatever the copy you downloaded was set to for versioning.

# Getting Started

Refer to the [vignette on installation](installation.html)
for details on installation and how to compile from source. An additional vignette
(TBD) [on optimization](optimization.html) is also included, discussing setting up a 
faster BLAS environment. Two ... blah blah.

# Contributing?

Examine our [contribution guidelines](Contributing.md) for more.
