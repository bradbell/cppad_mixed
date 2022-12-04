// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin ldlt_cholmod_sim_cov}
{xrst_spell
   cov
   covariance
   ll
   sim
}

Simulations with Covariance Corresponding to Factored Matrix
############################################################

Syntax
******
*ok* = *ldlt_obj* . ``sim_cov`` ( *w* , *v* )

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE
   // END_PROTOTYPE
}

Private
*******
The :ref:`ldlt_cholmod-name` class is an
:ref:`implementation detail<ldlt_cholmod@Private>` and not part of the
CppAD Mixed user API.

Purpose
*******
This function simulates a normal random vector with mean zero
and covariance :math:`H^{-1}` where

.. math::

   L D L^\R{T} = P H P^\R{T}

is the current factorization; see
:ref:`ldlt_cholmod@Factorization@H` ,
:ref:`ldlt_cholmod@Factorization@L` ,
:ref:`ldlt_cholmod@Factorization@D` , and
:ref:`ldlt_cholmod@Factorization@P` .

ldlt_obj
********
This object has prototype

   ``const CppAD::mixed::ldlt_cholmod`` *ldlt_obj*

In addition, it must have a previous call to
:ref:`ldlt_cholmod_update-name` .

w
*
This argument's
size is equal to the number of rows in :math:`H`.

v
*
This argument's
size is equal to the number of rows in :math:`H`.
The input value of its elements does not matter.
Upon return

.. math::

   v = P^\R{T} L^{-\R{T}} \tilde{D}^{-1/2} w

If :math:`w` is mean zero, variance identity white noise,
:math:`w \sim \B{N} ( 0 , I )`,
then :math:`v` will be mean zero and variance :math:`H^{-1}`,
:math:`v \sim \B{N} ( 0 , H^{-1} )`; see
:ref:`theory@Sparse Observed Information` .

Positive Definite
*****************
In the formula for :math:`v` above,
the matrix :math:`\tilde{D}` is a positive version of :math:`D`.
To be specific,

.. math::

   \tilde{D}_{i,i} = \left\{ \begin{array}{ll}
      D_{i,i} & \R{if} \; D_{i,i} \geq  \varepsilon^2 \; \max(D) \\
      \varepsilon^2 \; \max(D) & \R{otherwise}
   \end{array} \right.

where :math:`\varepsilon`
``std::numeric_limits<double>::epsilon()`` ,
and :math:`\max(D)` is the largest element in :math:`D`.

ok
**
If :math:`\max(D) > 0`, this routine terminates with *ok*
equal to true.
Otherwise it is false and the output values in *v*
are the same as their input values.

Example
*******
The file :ref:`ldlt_cholmod.cpp<ldlt_cholmod.cpp@sim_cov>` contains an
example and test that uses this function.

{xrst_end ldlt_cholmod_sim_cov}
*/
# include <cmath>
# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cppad/utility/index_sort.hpp>
# include <cassert>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN_PROTOTYPE
bool ldlt_cholmod::sim_cov(
   const CppAD::vector<double>& w  ,
   CppAD::vector<double>&       v  )
// END_PROTOTYPE
{  assert( update_called_ );
   //
   assert( factor_  != CPPAD_NULL );
   assert( rhs_     != CPPAD_NULL );
   assert( rhs_set_ != CPPAD_NULL );
   assert( rhs_->nrow == nrow_ );
   assert( rhs_->ncol == 1     );
   //
   assert( v.size() == nrow_ );
   //
   int*    factor_p  = (int *) factor_->p;
   double* factor_x  = (double *) factor_->x;
   int*    rhs_set_p = (int *) rhs_set_->p;
   int*    rhs_set_i = (int *) rhs_set_->i;
   double* rhs_x     = (double *) rhs_->x;
# ifndef NDEBUG
   int*    factor_i  = (int *) factor_->i;
# endif
   //
   // we will solve for all components on the right hand side
   rhs_set_p[0] = 0;
   rhs_set_p[1] = int( nrow_ );
   for(size_t i = 0; i < nrow_; i++)
      rhs_set_i[i] = int( i );
   //
   // set rhs_ = w
   for(size_t i = 0; i < nrow_; i++)
      rhs_x[i] = w[i];
   //
   // determine largest element in D
   double max_D = 0.0;
   for(size_t i = 0; i < nrow_; i++)
   {  // first element of each column is always the diagonal element
      assert( size_t( factor_i[ factor_p[i] ] ) == i );
      double di = factor_x[ factor_p[i] ];
      max_D     = std::max(max_D, di);
   }
   if( max_D <= 0.0 )
      return false;
   //
   // set rhs_ = \tilde{D}^{-1/2} w
   double eps = std::numeric_limits<double>::epsilon();
   eps        = eps * eps * max_D;
   for(size_t i = 0; i < nrow_; i++)
   {  // first element of each column is always the diagonal element
      assert( size_t( factor_i[ factor_p[i] ] ) == i );
      double di = factor_x[ factor_p[i] ];
      di        = std::max(di, eps);
      rhs_x[i]   = rhs_x[i] / std::sqrt( di );
   }
   // set sol_ = L^{-T} * D^{-1/2} w
   int sys = CHOLMOD_Lt;
# ifndef NDEBUG
   int flag =
# endif
   cholmod_solve2(
      sys,
      factor_,
      rhs_,
      rhs_set_,
      &sol_,
      &sol_set_,
      &work_one_,
      &work_two_,
      &common_
   );
   // check assumptions
   assert( flag == CHOLMOD_TRUE );
   assert( sol_set_->nrow == nrow_ );
   assert( sol_set_->ncol == 1 );
   assert( sol_set_->xtype == CHOLMOD_PATTERN );
   assert( sol_set_->packed == CHOLMOD_TRUE);
# ifndef NDEBUG
   int*   sol_set_p = (int *) sol_set_->p;
   assert( size_t( sol_set_p[1] ) == nrow_ );
# endif
   assert( sol_ != CPPAD_NULL   );
   assert( sol_->nrow == nrow_ );
   assert( sol_->ncol == 1     );
   //
   // set rhs_ = L^{-T} * D^{-1/2} w
   double* sol_x = (double *) sol_->x;
   for(size_t i = 0; i < nrow_; i++)
      rhs_x[i] = sol_x[i];
   //
   // set sol_ = P^T L^{-T} * D^{-1/2} w
   sys = CHOLMOD_Pt;
# ifndef NDEBUG
   flag =
# endif
   cholmod_solve2(
      sys,
      factor_,
      rhs_,
      rhs_set_,
      &sol_,
      &sol_set_,
      &work_one_,
      &work_two_,
      &common_
   );
   // check assumptions
   assert( flag == CHOLMOD_TRUE );
   assert( sol_set_->nrow == nrow_ );
   assert( sol_set_->ncol == 1 );
   assert( sol_set_->xtype == CHOLMOD_PATTERN );
   assert( sol_set_->packed == CHOLMOD_TRUE);
# ifndef NDEBUG
   sol_set_p = (int *) sol_set_->p;
   assert( size_t( sol_set_p[1] ) == nrow_ );
# endif
   assert( sol_ != CPPAD_NULL   );
   assert( sol_->nrow == nrow_ );
   assert( sol_->ncol == 1     );
   //
   // return v
   sol_x = (double *) sol_->x;
   for(size_t i = 0 ; i < nrow_; i++)
      v[i] = sol_x[i];
   //
   return true;
}

} } // END_CPPAD_MIXED_NAMESPACE
