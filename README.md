# Binomial Latent Subspace Estimation

Binomial LSE performs fast latent subspace estimation (LSE) for binomial data (genotypes) by performing fast principal component analysis (PCA) of single nucleotide
polymorphism (SNP) data. This implementation is based on the source code of [FlashPCA](https://github.com/gabraham/flashpca) and utilizes the [Spectra](https://github.com/yixuan/spectra/) and [Eigen](http://eigen.tuxfamily.org/) C++ libraries.

Latent subspace estimation described in [Chen and Storey 2015](https://arxiv.org/abs/1510.03497) is a modification of PCA that accounts for heteroskedasticity. We implement a scalable, low memory implementation of LSE for binomial data. This specific implementation uses the iteratively restarted Arnoldi method implemented by [FlashPCA](https://github.com/gabraham/flashpca) using [Spectra](https://github.com/yixuan/spectra/) to perform SVD/PCA.

## License
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This code is based on contributions from the following sources:
* [SparSNP](https://github.com/gabraham/SparSNP), Copyright (C) 2011-2012 Gad Abraham.
* [FlashPCA](https://github.com/gabraham/flashpca), Copyright (C) 2014-2020 Gad Abraham.
* [National ICT Australia](http://www.nicta.com.au).

## Building from source

To get the latest version:
   ```bash
   git clone git://github.com/alecmchiu/binomial_lse
   ```

### Requirements

On Linux:

* 64-bit OS
* g++ compiler
* [Eigen](http://eigen.tuxfamily.org), v3.2 or higher
   (if you get a compile error ``error: no match for 'operator/' in '1 / ((Eigen::MatrixBase...`` you'll need a more recent Eigen)
* [Spectra](https://github.com/yixuan/spectra/)
* [Boost](http://www.boost.org/), specifically boost_program_options/boost_program_options-mt.
* libgomp for openmp support

On Mac:

* [Homebrew](http://brew.sh) to install boost
* Eigen, as above
* Spectra, as above
* clang C++ compiler

### To install

The [Makefile](Makefile) contains three variables that need to be set according to where you have installed the Eigen
headers and Boost headers and libraries on your system. The default values for these are: 
   ```bash
   EIGEN_INC=/usr/local/include/eigen
   BOOST_INC=/usr/local/include/boost
   BOOST_LIB=/usr/local/lib
   SPECTRA_INC=spectra
   ```
   
 If your system has these libraries and header files in those locations, you can simply run make:
   ```bash
   cd binomial_lse
   make all
   ```
   
 If not, you can override their values on the make command line. For example,
 if you have the Eigen source in `/opt/eigen-3.2.5`, spectra headers in
 `/opt/spectra`, and Boost 1.59.0 installed into `/opt/boost-1.59.0`, you could run: 
   ```bash
   cd binomial_lse
   make all EIGEN_INC=/opt/eigen-3.2.5 \
      BOOST_INC=/opt/boost-1.59.0/include \
      BOOST_LIB=/opt/boost-1.59.0/lib \
      SPECTRA_INC=/opt/spectra
   ```

## Output

By default, binomial LSE produces the following files:

* `eigenvectors.txt`: the top k eigenvectors of the covariance matrix after adjustment for heteroskedasticity X X<sup>T</sup> / p - D. This is the file containing the subspace of interest.
* `pcs.txt`: the top k principal components (the projection of the data on the
eigenvectors, scaled by the eigenvalues
* `eigenvalues.txt`: the top k eigenvalues of X X<sup>T</sup> / p - D.
* `pve.txt`: the proportion of total variance explained by *each of the top k*
   eigenvectors (the total variance is given by the trace of the covariance
   matrix X X<sup>T</sup> / p - D, which is the same as the sum of all eigenvalues).
