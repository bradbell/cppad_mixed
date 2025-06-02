/*
// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-25 Bradley M. Bell
// ----------------------------------------------------------------------------
{xrst_begin sample_fixed}
{xrst_spell
  msg
  rng
  rcond
}

Sample Posterior for Fixed Effects
##################################

Syntax
******

| *error_msg* = *mixed_object* . ``sample_fixed`` (
| |tab| *sample* ,
| |tab| *hes_fixed_obj_rcv* ,
| |tab| *solution* ,
| |tab| *fixed_lower* ,
| |tab| *fixed_upper* ,
| |tab| *rcond*
| )

See Also
********
:ref:`sample_random-name`

Prototype
*********
{xrst_literal
   // BEGIN PROTOTYPE
   // END PROTOTYPE
}

Purpose
*******
This routine draws independent samples from
the asymptotic posterior distribution for the
fixed effects (given the model and the data).

Constant Fixed Effects
**********************
If the upper and lower limits for the *j*-th fixed effect are equal::

   fixed_lower[j] == fixed_upper[j]

we refer to the *j*-th fixed effect as constant.

Constraints
***********
Only the constant fixed effect constants are taken into account.
(This is an over estimate of the variance, but is faster to calculate
than trying to identify active constraints and treating them as equality
constraints.)
It can result in samples that are outside the limits
for the fixed effects that are not constant.
You may want to adjust the samples to be within
their upper and lower limits before you use them.

Covariance
**********
Each sample of the fixed effects
(excluding the constant fixed effects)
has covariance equal to the inverse of the
:ref:`information matrix<information_mat-name>`
(where the constant fixed effects have been removed).

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
We define::

   n_sample = sample_size / n_fixed

If *error_msg* is empty, upon return
for *i* = 0 , ..., *n_sample*-1 ,
*j* = 0 , ..., *n_fixed*-1 ::

   sample[i * n_fixed + j]

is the *j*-th component of the *i*-th sample of the
optimal fixed effects :math:`\hat{\theta}`.

#. These samples are independent for different :math:`i` .
#. For each :math:`i`, the corresponding sample has the specified
   :ref:`sample_fixed@Covariance` .
#. The constraints are only enforced for the
   :ref:`sample_fixed@Constant Fixed Effects` ; i.e.,
   if the *j*-th fixed effect is not constant,
   the sample may not satisfy the condition::

      fixed_lower[j] <= sample[i * n_fixed + j] <= fixed_upper[j]


hes_fixed_obj_rcv
*****************
This is a sparse matrix representation for the
lower triangle of the observed information matrix corresponding to
*solution* ; i.e., the matrix returned by

| *hes_fixed_obj_rcv* = *mixed_object* . ``hes_fixed_obj`` (
| |tab| *solution* , *random_opt*
| )

where *random_opt* is the optimal random effects corresponding
to *solution* .

solution
********
is the :ref:`optimize_fixed@solution`
for a the call to :ref:`optimize_fixed-name` corresponding to
*hes_fixed_obj_rcv* .
The only necessary information in this structure is
*solution* . ``fixed_opt`` .

fixed_lower
***********
is the same as
:ref:`optimize_fixed@fixed_lower`
in the call to ``optimize_fixed`` that corresponding to *solution* .

fixed_upper
***********
is the same as
:ref:`optimize_fixed@fixed_upper`
in the call to ``optimize_fixed`` that corresponding to *solution* .

Factorization
*************
The subset of the matrix *hes_fixed_obj_rcv* that corresponds to
fixed effects that are not constant is factored into L * D * L^T
where L is lower triangular and D is diagonal.

rcond
*****
This argument is optional and its input value does not matter.
Upon return it is the reciprocal of the condition number for
the diagonal matrix D in the factorization.
In other words, it is the minimum absolute entry in *D* divided
by the maximum absolute entry in D .
If the matrix D is singular, or any entry in D nan or infinite,
*rcond* is zero.

2DO: This result is not yet passing its test; see
:ref:`sample_fixed.cpp-name` .

error_msg
*********
If *error_msg* is empty (non-empty),
:ref:`sample_fixed@sample`
values have been calculated (have not been calculated).
If *error_msg* is non-empty,
it is a message describing the problem.
{xrst_toc_hidden
   example/user/sample_fixed.cpp
   src/eigen/sample_conditional.cpp
}
Example
*******
The file :ref:`sample_fixed.cpp-name` is an example
and test of ``sample_fixed`` .

Other Method
************
The routine :ref:`sample_conditional-name` , is a old method that is no
longer used for computing these samples.

{xrst_end sample_fixed}
------------------------------------------------------------------------------
*/
# include <Eigen/Core>
# include <Eigen/LU>
# include <Eigen/Cholesky>
# include <gsl/gsl_randist.h>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/manage_gsl_rng.hpp>
# include <cppad/mixed/undetermined.hpp>
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/configure.hpp>

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
// -------------------------------------------------------------------------
std::string cppad_mixed::try_sample_fixed(
   CppAD::vector<double>&                 sample               ,
   const d_sparse_rcv&                    hes_fixed_obj_rcv    ,
   const CppAD::mixed::fixed_solution&    solution             ,
   const CppAD::vector<double>&           fixed_lower          ,
   const CppAD::vector<double>&           fixed_upper          ,
   double&                                rcond                )
{
   // sample
   assert( sample.size() > 0 );
   assert( sample.size() % n_fixed_ == 0 );
   //
   // solution
   assert( solution.fixed_opt.size() == n_fixed_ );
   //
   // number of samples
   size_t n_sample = sample.size() / n_fixed_;
   //
   // optimal fixed effects
   const d_vector& fixed_opt( solution.fixed_opt );
   //
   // Determine the subset of variables that do not have lower equal to upper
   CppAD::vector<size_t> fixed2subset(n_fixed_);
   size_t n_subset = 0;
   for(size_t j = 0; j < n_fixed_; j++)
   {  if( fixed_lower[j] == fixed_upper[j] )
         fixed2subset[j] = n_fixed_;
      else
         fixed2subset[j] = n_subset++;
   }
   assert( n_subset <= n_fixed_ );
   //
   // sample
   if( n_subset == 0 )
   {  // All of the fixed effects have lower bound equal to upper bound
      for(size_t j = 0; j < n_fixed_; j++)
         for(size_t i_sample = 0; i_sample < n_sample; ++i_sample)
            sample[ i_sample * n_fixed_ + j] = fixed_lower[j];
      std::string error_msg = "";
      return error_msg;
   }
   //
   // create a sparse_rcv representation of information matrix on subset
   // and in column major order
   CppAD::vector<size_t> col_major = hes_fixed_obj_rcv.col_major();
   size_t count = 0;
   for(size_t ell = 0; ell < hes_fixed_obj_rcv.nnz(); ++ell)
   {  size_t k = col_major[ell];
      size_t i = hes_fixed_obj_rcv.row()[k];
      size_t j = hes_fixed_obj_rcv.col()[k];
      i        = fixed2subset[i];
      j        = fixed2subset[j];
      if( i != n_fixed_ && j != n_fixed_ && j <= i )
      {  assert( i < n_subset && j < n_subset );
         ++count;
      }
   }
   sparse_rc info_mat_rc(n_subset, n_subset, count);
   count = 0;
   for(size_t ell = 0; ell < hes_fixed_obj_rcv.nnz(); ++ell)
   {  size_t k = col_major[ell];
      size_t i = hes_fixed_obj_rcv.row()[k];
      size_t j = hes_fixed_obj_rcv.col()[k];
      i        = fixed2subset[i];
      j        = fixed2subset[j];
      if( i != n_fixed_ && j != n_fixed_ && j <= i )
      {  assert( i < n_subset && j < n_subset );
         info_mat_rc.set(count++, i, j);
      }
   }
   assert( count == info_mat_rc.nnz() );
   d_sparse_rcv info_mat_rcv( info_mat_rc );
   count = 0;
   for(size_t ell = 0; ell < hes_fixed_obj_rcv.nnz(); ++ell)
   {  size_t k = col_major[ell];
      size_t i = hes_fixed_obj_rcv.row()[k];
      size_t j = hes_fixed_obj_rcv.col()[k];
      double v = hes_fixed_obj_rcv.val()[k];
      i        = fixed2subset[i];
      j        = fixed2subset[j];
      if( i != n_fixed_ && j != n_fixed_ && j <= i )
      {  assert( i < n_subset && j < n_subset );
         info_mat_rcv.set(count++, v);
      }
   }
   assert( count == info_mat_rcv.nnz() );
   //
   // LDLT factorization of info_mat
   CPPAD_MIXED_LDLT_CLASS ldlt_info_mat(n_subset);
   ldlt_info_mat.init( info_mat_rcv.pat() );
   bool ok = ldlt_info_mat.update( info_mat_rcv );
   if( ! ok )
   {  std::string error_msg =
         "sample_fixed: fixed effects information matrix is singular";
      return error_msg;
   }
   size_t negative;
   ldlt_info_mat.logdet(negative);
   if( negative != 0 )
   {  std::string error_msg = "sample_fixed: "
         "fixed effects information matrix is not positive definite";
      return error_msg;
   }
   //
   // rcond
   rcond = ldlt_info_mat.rcond();
   //
   // -----------------------------------------------------------------------
   // Simulate the samples
   // -----------------------------------------------------------------------
   for(size_t i_sample = 0; i_sample < n_sample; i_sample++)
   {  d_vector w(n_subset);
      // simulate a normal with mean zero and variance one
      for(size_t k = 0; k < n_subset; k++)
         w[k] = gsl_ran_gaussian(CppAD::mixed::get_gsl_rng(), 1.0);
      //
      // set v to cholesky factor of info_mat^{-1} times w
      d_vector v(n_subset);
      ok = ldlt_info_mat.sim_cov(w, v);
      if( ! ok )
      {  std::string error_msg = "sample_fixed: fixed effects "
            " information matrix is not positive definite";
         return error_msg;
      }
      //
      // store in sample
      for(size_t j = 0; j < n_fixed_; j++)
      {  if( fixed2subset[j] == n_fixed_ )
            sample[ i_sample * n_fixed_ + j] = fixed_lower[j];
         else
         {  size_t k       = fixed2subset[j];
            //
            // store this component of the sample
            sample[ i_sample * n_fixed_ + j] = fixed_opt[j] + v[k];
         }
      }
   }
   // -----------------------------------------------------------------------
   return "";
}
// -------------------------------------------------------------------------
// BEGIN PROTOTYPE
std::string cppad_mixed::sample_fixed(
   CppAD::vector<double>&                 sample               ,
   const d_sparse_rcv&                    hes_fixed_obj_rcv    ,
   const CppAD::mixed::fixed_solution&    solution             ,
   const CppAD::vector<double>&           fixed_lower          ,
   const CppAD::vector<double>&           fixed_upper          ,
   double&                                rcond                )
// END PROTOTYPE
{  std::string error_msg = "";
   try
   {  error_msg = try_sample_fixed(
         sample            ,
         hes_fixed_obj_rcv ,
         solution          ,
         fixed_lower       ,
         fixed_upper       ,
         rcond
      );
   }
   catch(const std::exception& e)
   {  error_msg = "sample_fixed: std::exception: ";
      error_msg += e.what();
      return error_msg;
   }
   catch(const CppAD::mixed::exception& e)
   {  error_msg = e.message("sample_fixed");
      return error_msg;
   }
   return error_msg;
}
