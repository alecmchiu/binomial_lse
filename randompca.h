/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2014 Gad Abraham
 * All rights reserved.
 * 
 * FlashPCA code modified by Alec Chiu to implement binomial latent
 * subspace estimation as described by Chen and Storey 2015.
 */

#pragma once

#include <Eigen/QR>
#include <Eigen/SVD>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/fisher_f.hpp>

#include <Spectra/SymEigsSolver.h>

#include "data.h"

using namespace Eigen;

#define _8GiB 8589934592

#define METHOD_EIGEN 1
#define METHOD_SVD 2

#define MODE_PCA 1
#define MODE_CHECK_PCA 5

#define MEM_MODE_OFFLINE 1
#define MEM_MODE_ONLINE 2

#define LOWMEM 1
#define HIGHMEM 2

#define DIVISOR_NONE 0
#define DIVISOR_N1 1
#define DIVISOR_P 2

class RandomPCA {
   public:
      MatrixXd U, V, W, Px;
      VectorXd d;
      MatrixXd V0;
      double trace;
      VectorXd pve;
      MatrixXd X_meansd;
      ArrayXXd res;
      VectorXd err;
      double mse;
      double rmse;

      int stand_method_x;
      bool verbose;
      bool debug;
      int divisor;
      bool converged;

      void pca_fast(MatrixXd &X, unsigned int block_size,
	    unsigned int ndim,
	    unsigned int maxiter, double tol, long seed,
	    bool do_loadings);
      void pca_fast(Data &dat, unsigned int block_size,
	    unsigned int ndim,
	    unsigned int maxiter, double tol, long seed,
	    bool do_loadings);
      void check(Data& dat, unsigned int block_size,
	    std::string evec_file, std::string eval_file);
      void check(Data& dat, unsigned int block_size,
	    MatrixXd& evec, VectorXd& eval);
      void check(MatrixXd& X, MatrixXd& evec, VectorXd& eval);
      void project(Data& dat, unsigned int block_size);
};

