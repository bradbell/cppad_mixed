/*
// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
{xrst_begin hes_fixed_obj}

Compute the Hessian of The Fixed Effects Objective
##################################################

Syntax
******

| *hes_fixed_obj_rcv* = *mixed_object* . ``hes_fixed_obj`` (
| |tab| *fixed_vec* , *random_opt*
| )

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE
   // END_PROTOTYPE
}

Purpose
*******
Compute the Hessian of the fixed effects objective function
:math:`L ( \theta )`; see
:ref:`fixed effects objective<theory@Objective@Fixed Effects Objective, L(theta)>` .
The Hessian is

.. math::

   L^{(2)} ( \hat{\theta} )

Absolute value terms in the
:ref:`problem@Negative Log-Density Vector`
for the :ref:`fix_likelihood-name` are not include in this Hessian
(because they do not have a derivative, let alone Hessian, at zero).

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

fixed_vec
*********
is the vector of fixed effects :math:`\theta` at which
the Hessian is evaluated.

random_opt
**********
is the optimal random effects corresponding to this value for the
fixed effects; i.e.,
:ref:`u^(theta)<theory@Optimal Random Effects, u^(theta)>` .

hes_fixed_obj_rcv
*****************
The return value is a
:ref:`typedef@Sparse Types@d_sparse_rcv` representation
of the lower triangle of the Hessian.
(The Hessian is symmetric and hence determined by its lower triangle.)
Absolute value terms in the
:ref:`problem@Negative Log-Density Vector`
for the :ref:`fix_likelihood-name` are not include in this Hessian
because they do not have a derivative (let alone Hessian) at zero.
{xrst_toc_hidden
   example/user/hes_fixed_obj.cpp
}

Example
*******
The file :ref:`hes_fixed_obj.cpp-name` contains an example and
test of this routine. It returns true for success and false for failure.

{xrst_end hes_fixed_obj}
------------------------------------------------------------------------------
*/
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/triple2eigen.hpp>
# include <cppad/mixed/exception.hpp>

CppAD::mixed::d_sparse_rcv cppad_mixed::try_hes_fixed_obj(
   const d_vector& fixed_vec             ,
   const d_vector& random_opt            )
{
   typedef Eigen::SparseMatrix<double, Eigen::ColMajor>      eigen_sparse;
   typedef eigen_sparse::InnerIterator                       sparse_itr;
   //
   assert( fixed_vec.size()  == n_fixed_ );
   assert( random_opt.size() == n_random_ );
   //
   // ----------------------------------------------------------------------
   // Hessian w.r.t. fixed effect of laplace part of fixed evects objective
   eigen_sparse laplace_hes;
   if( n_random_ == 0 )
      laplace_hes.resize( int(n_fixed_), int(n_fixed_) );
   else
   {  // ------------------------------------------------------------------
      // update the cholesky factor for this fixed and random effect
      update_factor(fixed_vec, random_opt);
      // ------------------------------------------------------------------
      // If Quasi-Newton method was used, must initilaize routines
      // that are only used for the Hessian calculation; see initilaize.cpp
      if( ! init_laplace_obj_fun_done_ )
      {  assert( ! init_laplace_obj_done_ );
         init_laplace_obj(fixed_vec, random_opt);
         assert( init_laplace_obj_done_ );
      }
      // ------------------------------------------------------------------
      // Lower triangle of Hessian w.r.t. fixed effects
      // for laplace part of the fixed effect objective, no constraints
      d_vector w_laplace(A_rcv_.nr() + 1);
      w_laplace[0] = 1.0;
      for(size_t j = 1; j <=A_rcv_.nr(); j++)
         w_laplace[j] = 0.0;
      //
      CppAD::mixed::sparse_mat_info laplace_info;
      laplace_obj_hes(
         fixed_vec,
         random_opt,
         w_laplace,
         laplace_info.row,
         laplace_info.col,
         laplace_info.val
      );
      CppAD::mixed::triple2eigen(
         laplace_hes       ,
         n_fixed_          ,
         n_fixed_          ,
         laplace_info.row  ,
         laplace_info.col  ,
         laplace_info.val
      );
   }
   // --------------------------------------------------------------------
   // Lower triangle of Hessian of the fixed part of fixed effects objective
   eigen_sparse fixed_only_hes;
   size_t n_fix_like = 0;
   if( fix_like_fun_.size_var() == 0 )
      fixed_only_hes.resize( int(n_fixed_) , int(n_fixed_) );
   else
   {  n_fix_like = fix_like_fun_.Range();
      d_vector w_fix(n_fix_like);
      w_fix[0] = 1.0;
      for(size_t j = 1; j < n_fix_like; j++)
         w_fix[0] = 0.0;
      //
      CppAD::mixed::sparse_mat_info fix_info;
      fix_like_hes(
            fixed_vec, w_fix, fix_info.row, fix_info.col, fix_info.val
      );
      CppAD::mixed::triple2eigen(
         fixed_only_hes  ,
         n_fixed_        ,
         n_fixed_        ,
         fix_info.row    ,
         fix_info.col    ,
         fix_info.val
      );
   }
   // -----------------------------------------------------------------------
   // Hessian of fixed effects objective (observed information matrix)
   eigen_sparse total_hes = laplace_hes + fixed_only_hes;
   //
   // convert from eigen sparse matrix to d_sparse_rcv
   size_t nr  = total_hes.rows();
   size_t nc  = total_hes.cols();
   size_t nnz = total_hes.nonZeros();
   sparse_rc pattern(nr, nc, nnz);
   d_vector val(nnz);
   size_t k = 0;
   for(size_t j = 0; j < n_fixed_; j++)
   {  for(sparse_itr itr(total_hes, int(j) ); itr; ++itr)
      {  assert( size_t( itr.col() ) == j );
         size_t i = itr.row();
         pattern.set(k, i, j);
         val[k] = itr.value();
         ++k;
      }
   }
   assert( k == nnz );
   d_sparse_rcv total_rcv(pattern);
   for(k = 0;  k < nnz; k++)
      total_rcv.set(k, val[k]);
   //
   return total_rcv;
}
// BEGIN_PROTOTYPE
CppAD::mixed::d_sparse_rcv cppad_mixed::hes_fixed_obj(
   const d_vector& fixed_vec            ,
   const d_vector& random_opt           )
// END_PROTOTYPE
{  d_sparse_rcv result;
   try
   {  result = try_hes_fixed_obj(fixed_vec, random_opt);
   }
   catch(const std::exception& e)
   {  std::string error_message = "hes_fixed_obj: std::exception: ";
      error_message += e.what();
      fatal_error(error_message);
      assert(false);
   }
   catch(const CppAD::mixed::exception& e)
   {  std::string error_message = e.message("hes_fixed_obj");
      fatal_error(error_message);
      assert(false);
   }
   //
   return result;
}
