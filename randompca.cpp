/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Copyright (C) 2014-2016 Gad Abraham
 * All rights reserved.
 * 
 * FlashPCA code modified by Alec Chiu to implement binomial latent
 * subspace estimation as described by Chen and Storey 2015.
 */

#include "randompca.h"
#include "util.h"
#include "svdwide.h"

void RandomPCA::pca_fast(MatrixXd& X, unsigned int block_size,
   unsigned int ndim, unsigned int maxiter,
   double tol, long seed, bool do_loadings)
{
   unsigned int N, p;

   X_meansd = standardise(X, stand_method_x, verbose);
   N = X.rows();
   p = X.cols();

   SVDWide op(X, verbose);
   Spectra::SymEigsSolver<double,
      Spectra::LARGEST_ALGE, SVDWide> eigs(&op, ndim, ndim * 2 + 1);

   eigs.init();
   eigs.compute(maxiter, tol);

   double div = 1;
   if(divisor == DIVISOR_N1)
      div = N - 1;
   else if(divisor == DIVISOR_P)
      div = p;

   if(eigs.info() == Spectra::SUCCESSFUL)
   {
      U = eigs.eigenvectors();
      // Note: _eigenvalues_, not singular values
      d = eigs.eigenvalues().array() / div;
      if(do_loadings)
      {
         VectorXd s = d.array().sqrt().inverse() / sqrt(div);
         V.noalias() = X.transpose() * U * s.asDiagonal();
      }
      trace = X.array().square().sum() / div;
      pve = d / trace;
      Px = U * d.array().sqrt().matrix().asDiagonal();

      verbose && STDOUT << timestamp() << "GRM trace: " << trace << std::endl;
   }
   else
   {
      throw new std::runtime_error(
	 std::string("Spectra eigen-decomposition was not successful")
	    + ", status: " + std::to_string(eigs.info()));
   }
}

void RandomPCA::pca_fast(Data& dat, unsigned int block_size,
   unsigned int ndim, unsigned int maxiter, double tol,
   long seed, bool do_loadings)
{
   unsigned int N = dat.N, p = dat.nsnps;
   SVDWideOnline op(dat, block_size, stand_method_x, verbose);
   Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE,
      SVDWideOnline> eigs(&op, ndim, ndim * 2 + 1);

   eigs.init();
   eigs.compute(maxiter, tol);

   double div = 1;
   if(divisor == DIVISOR_N1)
      div = N - 1;
   else if(divisor == DIVISOR_P)
      div = p;

   if(eigs.info() == Spectra::SUCCESSFUL)
   {
      U = eigs.eigenvectors();
      // Note: _eigenvalues_, not singular values
      d = eigs.eigenvalues().array() / div;
      if(do_loadings)
      {
         V = MatrixXd::Zero(dat.nsnps, U.cols());
         verbose && STDOUT << "Computing loadings" << std::endl;
         VectorXd v(dat.nsnps);
         for(unsigned int j = 0 ; j < U.cols() ; j++)
         {
	    verbose && STDOUT << "loading " << j << std::endl;
            VectorXd u = U.col(j);
            op.crossprod(u.data(), v.data());
            double s = d(j);
            V.col(j) = v * (1.0 / sqrt(s)) / sqrt(div);
         }
      }
      trace = op.trace / div;
      pve = d / trace;
      Px = U * d.array().sqrt().matrix().asDiagonal();
      X_meansd = dat.X_meansd; // TODO: duplication

      verbose && STDOUT << timestamp() << "GRM trace: " << trace << std::endl;
   }
   else
   {
      throw new std::runtime_error(
	 std::string("Spectra eigen-decomposition was not successful")
	    + ", status: " + std::to_string(eigs.info()));
   }
}

void RandomPCA::check(Data& dat, unsigned int block_size,
   std::string evec_file, std::string eval_file)
{
   // Read eigenvalues
   // Expects no header, no rownames, one eigenvalue per row
   verbose && STDOUT << timestamp() << "Loading eigenvalue file '"
       << eval_file << "'" << std::endl;
   NamedMatrixWrapper M1 = read_text(eval_file.c_str(), 1, -1, 0);
   MatrixXd ev = M1.X;
   if(ev.rows() == 0)
      throw std::runtime_error("No eigenvalues found in file");
   VectorXd eval = ev.col(0);

   // Read eigenvectors
   // Expects header (colnames), FID and IID cols
   verbose && STDOUT << timestamp() << "Loading eigenvector file '"
       << evec_file << "'" << std::endl;
   NamedMatrixWrapper M2 = read_text(evec_file.c_str(), 3, -1, 1);
   MatrixXd& evec = M2.X;

   if(evec.rows() != dat.N)
      throw std::runtime_error(
	 std::string("Eigenvector dimension doesn't match data dimension")
	    + " (evec.rows = " + std::to_string(evec.rows())
	    + "; dat.N = " + std::to_string(dat.N) + ")");

   if(eval.size() != evec.cols())
      throw std::runtime_error(
	 "Eigenvector dimension doesn't match the number of eigenvalues");

   check(dat, block_size, evec, eval);
}

// Check the eigenvalues/eigenvectors, computing the root mean squared error
// of E = X X' U / div - U D^2, i.e., averaged over the n * k dimensions.
// assumes WIDE
void RandomPCA::check(Data& dat, unsigned int block_size,
   MatrixXd& evec, VectorXd& eval)
{
   SVDWideOnline op(dat, block_size, 1, verbose);

   unsigned int K = std::min(evec.cols(), eval.size());

   // X X' U / div = U D^2
   verbose && STDOUT << timestamp()
      << "Checking mean square error between (X X' U) / div and (U D^2)"
      << " for " << K << " dimensions"
      << std::endl;

   double div = 1;
   if(divisor == DIVISOR_N1)
      div = dat.N - 1;
   else if(divisor == DIVISOR_P)
      div = dat.nsnps;

   MatrixXd XXU = op.perform_op_mat(evec);
   XXU /= div;
   MatrixXd UD2 = evec * eval.asDiagonal();

   RowVectorXd rerr = (XXU - UD2).colwise().squaredNorm();
   err = rerr.transpose();

   for(unsigned int j = 0 ; j < K ; j++)
   {
      verbose && STDOUT << timestamp() << "eval(" << (j + 1)
         << "): " << eval(j) << ", sum squared error: "
         << err(j) << std::endl;
   }

   mse = err.sum() / (dat.N * K);
   rmse = std::sqrt(mse);

   verbose && STDOUT << timestamp() << "Mean squared error: " << mse
      << ", Root mean squared error: " << rmse
      << " (n=" << dat.N << ")" << std::endl;

}

void RandomPCA::check(MatrixXd& X, MatrixXd& evec, VectorXd& eval)
{
   X_meansd = standardise(X, stand_method_x, verbose);
   unsigned int K = std::min(evec.cols(), eval.size());
   unsigned int N = X.rows();
   unsigned int p = X.cols();

   // X X' U / div = U D^2
   verbose && STDOUT << timestamp()
      << "Checking mean square error between (X X' U) and (U D^2)"
      << " for " << K << " dimensions"
      << std::endl;

   double div = 1;
   if(divisor == DIVISOR_N1)
      div = N - 1;
   else if(divisor == DIVISOR_P)
      div = p;

   MatrixXd XXU = X * (X.transpose() * evec) / div;
   MatrixXd UD2 = evec * eval.asDiagonal();

   RowVectorXd rerr = (XXU - UD2).colwise().squaredNorm();
   err = rerr.transpose();

   for(unsigned int j = 0 ; j < K ; j++)
   {
      verbose && STDOUT << timestamp() << "eval(" << (j + 1)
         << "): " << eval(j) << ", sum squared error: "
         << err(j) << std::endl;
   }

   mse = err.sum() / (N * K);
   rmse = std::sqrt(mse);

   verbose && STDOUT << timestamp() << "Mean squared error: " << mse
      << ", Root mean squared error: " << rmse
      << " (n=" << N << ")" << std::endl;
}

// Project new samples onto existing principal components.
//
// Doesn't do a lot of sanity checking.
//
// Assumes:
// - The loadings matrix V must have been set already
// - dat.X_meansd has beeen set
// - dat.use_preloaded_maf has been set, if needed
void RandomPCA::project(Data& dat, unsigned int block_size)
{
   // Check that the SNPs in the data match

   SVDWideOnline op(dat, block_size, 1, verbose);

   unsigned int k = V.cols();
   Px = MatrixXd::Zero(dat.N, k);

   double div = 1;
   if(divisor == DIVISOR_N1)
      div = dat.N - 1;
   else if(divisor == DIVISOR_P)
      div = V.rows();

   VectorXd pxi(dat.N);
   for(unsigned int i = 0 ; i < k ; i++)
   {
      MatrixXd v = V.col(i);
      op.prod(v.data(), pxi.data());
      Px.col(i) = pxi.array() / sqrt(div); // X V = U D
   }
}

