/*
// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
{xrst_begin sample_conditional}
{xrst_spell
  rng
}

Sample Posterior for Fixed Effects Using Conditional Covariance
###############################################################

Syntax
******

| *mixed_object* . ``sample_conditional`` (
| |tab| *sample* ,
| |tab| *information_info* ,
| |tab| *solution* ,
| |tab| *fixed_lower* ,
| |tab| *fixed_upper* ,
| |tab| *random_opt*
| )

Replaced
********
This routine uses the
:ref:`sample_conditional@Theory@Conditional Covariance`
to sample the fixed effects.
This required inverting the :ref:`information_mat-name` .
It has been replaced by using the :ref:`sample_fixed-name` because it is faster.

Prototype
*********
{xrst_literal
   // BEGIN PROTOTYPE
   // END PROTOTYPE
}

Purpose
*******
This routine draw samples from
the asymptotic posterior distribution for the
optimal fixed effects (given the model and the data).

manage_gsl_rng
**************
It is assumed that
:ref:`manage_gsl_rng@get_gsl_rng` will return
a pointer to a GSL random number generator.

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

sample
******
The size *sample* . ``size`` () is a multiple of
:ref:`derived_ctor@n_fixed` .
The input value of its elements does not matter.
We define

   *n_sample* = *sample_size* / *n_fixed*

Upon return,
for ``i`` = 0 , ..., ``n_sample`` *-1* ,
``j`` = 0 , ..., ``n_fixed`` *-1* ,

   *sample* [ *i* * *n_fixed* + *j*  ]

is the *j*-th component of the *i*-th sample of the
optimal fixed effects :math:`\hat{\theta}`.
These samples are independent for different :math:`i`,
and for fixed :math:`i`, they have the
:ref:`sample_conditional@Theory@Conditional Covariance`
:math:`D`.

information_info
****************
This is a sparse matrix representation for the
lower triangle of the observed information matrix corresponding to
*solution* ; i.e., the matrix returned by

| *information_info* = *mixed_object* . ``information_mat`` (
| |tab| *solution* , *random_options* , *random_lower* , *random_upper* , *random_in*
| )

solution
********
is the :ref:`optimize_fixed@solution`
for a the call to :ref:`optimize_fixed-name` corresponding to
*information_info* .

fixed_lower
***********
is the same as
:ref:`optimize_fixed@fixed_lower`
in the call to ``optimize_fixed`` that corresponds to *solution* .

fixed_upper
***********
is the same as
:ref:`optimize_fixed@fixed_upper`
in the call to ``optimize_fixed`` that corresponds to *solution* .

random_opt
**********
is the optimal random effects corresponding to the solution; i.e.

| |tab| *random_opt* = *mixed_object* . ``optimize_random`` (
| |tab| |tab| *random_options* ,
| |tab| |tab| *solution* . ``fixed_opt`` ,
| |tab| |tab| *random_lower* ,
| |tab| |tab| *random_upper* ,
| |tab| |tab| *random_in*
| |tab| )

*random_options* ,
*random_lower* ,
*random_upper* , and
*random_in* , are the same
as in the call to ``optimize_fixed`` that corresponds to *solution* .

Theory
******

Notation
========
Given two random vectors :math:`u` and :math:`v`,
we use the notation :math:`\B{C}( u , v )`
for the corresponding covariance matrix;
i.e.,

.. math::

   \B{C}( u , v )
   =
   \B{E} \left( [ u - \B{E} (u) ] [ v - \B{E} (v) ]^\R{T} \right)

Fixed Effects Subset
====================
We use :math:`\alpha` for the vector of fixed effects that do not have
their upper or lower bound active (or equal); i.e., if *j* is such that

   *solution* . ``fixed_lag`` [ *j* ] == 0.0 && *fixed_lower* [ *j* ] < *fixed_upper* [ *j* ]

then :math:`\theta_j` is one of the components in :math:`\alpha`.
Note that each value of :math:`\alpha` has a corresponding value for
:math:`\theta` where the active bounds are used for the components
not in :math:`\alpha`.

Unconstrained Subset Covariance
===============================
Note that the bound constraints do not apply to the subset of fixed effects
represented by :math:`\alpha`.
We use :math:`\tilde{L} ( \alpha )` to denote the
:ref:`fixed effects objective<theory@Objective@Fixed Effects Objective, L(theta)>`
as a function of :math:`\alpha` and
where the absolute values terms in :ref:`fix_likelihood-name` are excluded.
We use :math:`\tilde{\alpha}` for the unconstrained optimal estimate
of the subset of fixed effects and
approximate its auto-covariance by

.. math::

   \B{C} ( \tilde{\alpha} , \tilde{\alpha} )
   =
   H^{-1}

Here :math:`H` is the Hessian corresponding to *information_info* .
Note that *information_info* is the observed information matrix
corresponding to all the fixed effects :math:`\theta`.

Constraint Equations
====================
Let :math:`n` be the number of fixed effects in :math:`\alpha`,
:math:`m` the number of active constraints (not counting bounds),
and the equations :math:`e( \alpha ) = b` be those active constraints.
Here :math:`e : \B{R}^n \rightarrow \B{R}^m` and :math:`b \in \B{R}^m`
and the inequality constraints have been converted to equalities at the
active bounds (excluding the bounds on the fixed effects).
Define the random variable the approximation for :math:`e( \alpha )` by

.. math::

   \tilde{e} ( \alpha ) =
   e \left( \hat{\alpha} \right) + e^{(1)} \left( \hat{\alpha} \right)
      \left( \alpha - \hat{\alpha} \right)

where :math:`\hat{\alpha}` is the subset of the optional estimate
for the fixed effects *solution* . ``fixed_opt`` .

Conditional Covariance
======================
We approximate the distribution for
:math:`\tilde{\alpha}` normal,
and the distribution for :math:`\hat{\alpha}`
as the conditional distribution of :math:`\tilde{\alpha}` given
the value of :math:`\tilde{e} ( \tilde{\alpha} )`; i.e.,

.. math::

   \B{C} \left( \hat{\alpha} \W{,} \hat{\alpha} \right)
   =
   \B{C} \left( \tilde{\alpha} \W{,} \tilde{\alpha} \right)
   -
   \B{C} \left( \tilde{\alpha} \W{,} \tilde{e} \right)
   \B{C} \left( \tilde{e}  \W{,} \tilde{e} \right)^{-1}
   \B{C} \left( \tilde{e}  \W{,} \tilde{\alpha} \right)

Using the notation
:math:`D = \B{C} \left( \hat{\alpha} \W{,} \hat{\alpha} \right)`,
:math:`C = \B{C} \left( \tilde{\alpha} \W{,} \tilde{\alpha} \right)`,
:math:`E = e^{(1)} \left( \hat{\alpha} \right)`,
we have

.. math::

   D = C - C E^\R{T} \left( E C E^\R{T} \right)^{-1}  E C

Example
*******
The file :ref:`sample_fixed.cpp-name` is an example
and test of ``sample_conditional`` was used before it was
:ref:`sample_conditional@Replaced` .

{xrst_end sample_conditional}
------------------------------------------------------------------------------
*/
# include <Eigen/Core>
# include <Eigen/LU>
# include <Eigen/Cholesky>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/manage_gsl_rng.hpp>
# include <gsl/gsl_randist.h>

# define DEBUG_PRINT 0

namespace {
   using Eigen::Dynamic;
   using CppAD::mixed::get_gsl_rng;
   typedef Eigen::Matrix<double, Dynamic, Dynamic>     double_mat;
   typedef Eigen::Matrix<double, Dynamic, 1>           double_vec;
   typedef Eigen::Matrix<size_t, Dynamic, 1>           size_vec;
   typedef Eigen::LDLT<double_mat, Eigen::Lower>       double_cholesky;
   typedef Eigen::PermutationMatrix<Dynamic, Dynamic, int>  permutation_mat;
   //
# if DEBUG_PRINT
   void print(const char* name , const double_mat& mat)
   {  std::cout << "\n" << name << " =\n" << mat << "\n"; }
   void print(const char* name , double_vec& vec)
   {  std::cout << "\n" << name << "^T = " << vec.transpose() << "\n"; }
   void print(const char* name , size_vec& vec)
   {  std::cout << "\n" << name << "^T = " << vec.transpose() << "\n"; }
# endif
}

// BEGIN PROTOTYPE
void cppad_mixed::sample_conditional(
   CppAD::vector<double>&                 sample               ,
   const CppAD::mixed::sparse_mat_info&   information_info     ,
   const CppAD::mixed::fixed_solution&    solution             ,
   const CppAD::vector<double>&           fixed_lower          ,
   const CppAD::vector<double>&           fixed_upper          ,
   const CppAD::vector<double>&           random_opt           )
// END PROTOTYPE
{
   // number of fixed constraints
   size_t n_fix_con = 0;
   if( fix_con_fun_.size_var() != 0 )
      n_fix_con = fix_con_fun_.Range();
   //
   // sample
   assert( sample.size() > 0 );
   assert( sample.size() % n_fixed_ == 0 );
   // information_info
   assert( information_info.row.size() == information_info.col.size() );
   assert( information_info.row.size() == information_info.val.size() );
   // solution
   assert( solution.fixed_opt.size() == n_fixed_ );
   assert( solution.fixed_lag.size() == n_fixed_ );
   assert( solution.fix_con_lag.size() == n_fix_con );
   assert( solution.ran_con_lag.size() == A_rcv_.nr() );
   // random_opt
   assert( random_opt.size() == n_random_ );
   //
   // number of samples
   size_t n_sample = sample.size() / n_fixed_;
   //
   // optimal fixed effects
   const d_vector& fixed_opt( solution.fixed_opt );
   // -----------------------------------------------------------------------
   // update the cholesky factor for this fixed and random effect
   if( n_random_ > 0 )
      update_factor(fixed_opt, random_opt);
   // -----------------------------------------------------------------------
   // subset of variables that do not have active bounds
   // mapping from fixed index to subset index and back
   CppAD::vector<size_t> fixed2subset(n_fixed_);
   size_t n_subset = 0;
   for(size_t j = 0; j < n_fixed_; j++)
   {  if( solution.fixed_lag[j] != 0.0 )
         fixed2subset[j] = n_fixed_;
      else if( fixed_lower[j] == fixed_upper[j] )
         fixed2subset[j] = n_fixed_;
      else
         fixed2subset[j] = n_subset++;
   }
   assert( n_subset <= n_fixed_ );
   // -----------------------------------------------------------------------
   // Create con_mat
   //
   // number of variables in the subset that do not have active bounds
   //
   // number fixed constraints active
   size_t n_fix_active = 0;
   size_vec fix_active_index(n_fix_con);
   for(size_t i = 0; i < n_fix_con; i++)
   {  fix_active_index[i] = n_fixed_;
      if( solution.fix_con_lag[i] != 0.0 )
         fix_active_index[i] = n_fix_active++;
   }
   // number of random constraints active
   size_t n_ran_active = A_rcv_.nr();
   //
   // matrix with all the active constraints
   // (not counting bound constraints)
   size_t n_con_active = n_fix_active + n_ran_active;
   double_mat con_mat = double_mat::Zero(n_con_active, n_subset);
   size_t con_row = 0;
   //
   // put fixed constraints in con_mat
   if( n_fix_con > 0 )
   {  // jacobian of the fixed constraints
      CppAD::mixed::sparse_mat_info fix_con_info;
      fix_con_jac(
         fixed_opt, fix_con_info.row, fix_con_info.col, fix_con_info.val
      );
      size_t K = fix_con_info.row.size();
      for(size_t k = 0; k < K; k++)
      {  size_t r  = fix_con_info.row[k];
         size_t c  = fix_con_info.col[k];
         double v  = fix_con_info.val[k];
         //
         bool r_active  = solution.fix_con_lag[r] != 0.0;
         bool in_subset = fixed2subset[c] != n_fixed_;
         if( r_active & in_subset )
         {  assert( fix_active_index[r] != n_fixed_ );
            size_t i  = fix_active_index[r];
            size_t j  = fixed2subset[c];
            con_mat(con_row + i, j) = v;
         }
      }
      con_row += n_fix_active;
   }
   // put random constraints in con_mat
   CppAD::mixed::sparse_mat_info ran_con_info;
   if( n_ran_active > 0 )
   {  assert( con_row == n_fix_active );
      //
      // jacobian of the random constraints
      // sparsity pattern
      ran_con_jac(fixed_opt, random_opt, ran_con_info);
      // values
      ran_con_jac(fixed_opt, random_opt, ran_con_info);
      //
      size_t K = ran_con_info.row.size();
      for(size_t k = 0; k < K; k++)
      {  size_t r       = ran_con_info.row[k];
         size_t c       = ran_con_info.col[k];
         double v       = ran_con_info.val[k];
         bool in_subset = fixed2subset[c] != n_fixed_;
         if( in_subset )
         {  size_t j  = fixed2subset[c];
            con_mat(con_row + r, j) = v;
         }
      }
   }
   // -----------------------------------------------------------------------
   // compute conditional covariance
   // ----------------------------------------------------------------------
   //
   // information matrix
   double_mat info_mat = double_mat::Zero(n_subset, n_subset);
   for(size_t k = 0; k < information_info.row.size(); k++)
   {  // note only lower triangle is stored in information_info
      size_t r = information_info.row[k];
      size_t c = information_info.col[k];
      double v = information_info.val[k];
      //
      bool in_subset = fixed2subset[r] != n_fixed_;
      in_subset     &= fixed2subset[c] != n_fixed_;
      //
      if( in_subset )
      {  size_t i = fixed2subset[r];
         size_t j = fixed2subset[c];
         info_mat(i, j) = v;
         info_mat(j, i) = v;
      }
   }
   //
   // covariance matrix with out constraints
   double_mat C = double_mat( info_mat.inverse() );
   //
   // conditional covariance
   double_mat& E(con_mat);
   double_mat  EC    = E * C;
   double_mat  ECET  = EC * E.transpose();
   double_mat  D     = C - EC.transpose() * ECET.inverse() * EC;
   //
   // LDLT factorizaton of D
   double_cholesky cholesky;
   cholesky.compute(D);
   //
   // diagonal elements of LDLT factorization
   double_vec diag      = cholesky.vectorD();
   double_vec diag_root(n_subset);
   for(size_t j = 0; j < n_subset; j++)
   {  if( diag[j] > 0.0 )
         diag_root[j] = std::sqrt( diag[j] );
      else
         diag_root[j] = 0.0;
   }
   double_mat L      = cholesky.matrixL();
   permutation_mat P = permutation_mat( cholesky.transpositionsP() );
   // -----------------------------------------------------------------------
   // Simulate the samples
   // -----------------------------------------------------------------------
   for(size_t i_sample = 0; i_sample < n_sample; i_sample++)
   {  double_vec w(n_subset);
      // simulate a normal with mean zero and variance sqrt{D(k,k)}
      for(size_t k = 0; k < n_subset; k++)
         w[k] = diag_root[k] * gsl_ran_gaussian(get_gsl_rng(), 1.0);
      // multily by Cholesky factor
      double_vec s = P.transpose() * L * w;
      //
      // store this sample
      for(size_t j = 0; j < n_fixed_; j++)
      {  if( fixed2subset[j] == n_fixed_ )
            sample[ i_sample * n_fixed_ + j] = fixed_opt[j];
         else
         {  size_t k       = fixed2subset[j];
            double fixed_j = fixed_opt[j] + s[k];
            //
            // check if this component went out of bounds
            fixed_j = std::min(fixed_j, fixed_upper[j]);
            fixed_j = std::max(fixed_j, fixed_lower[j]);
            //
            // store this component of the sample
            sample[ i_sample * n_fixed_ + j] = fixed_j;
         }
      }
   }
   // -----------------------------------------------------------------------
   return;
}
