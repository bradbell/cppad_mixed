// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>
/*
{xrst_begin laplace_obj_hes dev}
{xrst_spell
  nr
}

Hessian of Laplace Objective and Random Constraints
###################################################

Syntax
******

| *mixed_object* . ``laplace_obj_hes`` (
| |tab| *fixed_vec* , *random_vec* , *weight* , *row_out* , *col_out* , *val_out*
| )

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

Purpose
*******
This routine computes the Hessian, w.r.t the fixed effects,
of a vector times the random part of the of the objective
and random constraints; i.e.,

.. math::

   w_0 H_{\beta, \beta} ( \beta , \theta , u )
   +
   \sum_{i=1}^m
   w_i B^{i-1}_{\beta, \beta} ( \beta , \theta , u )

where :math:`m` is the number of rows in the
:ref:`random constraint matrix<problem@Notation@Random Constraint Matrix, A>` ;
see
:ref:`H(beta, theta, u)<theory@Approximate Laplace Objective, H(beta, theta, u)>`
and
:ref:`B(beta, theta, u)<theory@Approximate Random Constraint Function, B(beta, theta, u)>` .

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

fixed_vec
*********
This argument has prototype

   ``const CppAD::vector<double>&`` *fixed_vec*

It specifies the value of the
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>`
vector :math:`\theta`.

random_vec
**********
This argument has prototype

   ``const CppAD::vector<double>&`` *random_vec*

Given the fixed effects :math:`\theta`, it is the corresponding
:ref:`optimal random effects<theory@Optimal Random Effects, u^(theta)>`
:math:`\hat{u} ( \theta )`.

weight
******
This argument has prototype

   ``const CppAD::vector<double>&`` *weight*

It determines the weighting vector :math:`w` in the summation

.. math::

   w_0 H_{\beta, \beta} ( \beta , \theta , u )
   +
   \sum_{i=1}^m
   w_i B^{i-1}_{\beta, \beta} ( \beta , \theta , u )

The size of *weight* is
:ref:`A_rcv.nr()<derived_ctor@A_rcv>` plus one.

row_out
*******
This argument has prototype

   ``CppAD::vector<size_t>&`` *row_out*

If the input size of this array is non-zero,
the entire vector must be the same
as for a previous call to ``laplace_obj_hes`` .
If it's input size is zero,
upon return it contains the row indices for the Hessian elements
that are possibly non-zero;

   *row_out* [ *k* ] < *n_fixed*

for all *k* = 0 , ..., *row_out* . ``size`` () ``-1``

col_out
*******
This argument has prototype

   ``CppAD::vector<size_t>&`` *col_out*

If the input size of this array is non-zero,
the entire vector must be the same as for
a previous call to ``laplace_obj_hes`` .
If it's input size is zero,
upon return it contains the column indices for the Hessian elements
that are possibly non-zero (and will have the same size as *row_out* ).
Note that only the lower triangle of the Hessian is computed and hence

   *col_out* [ *k* ] <= *row_out* [ *k* ]

for all *k* = 0 , ..., *row_out* . ``size`` () ``-1``

val_out
*******
This argument has prototype

   ``CppAD::vector<double>&`` *val_out*

If the input size of this array is non-zero, it must have the same size
as for a previous call to ``laplace_obj_hes`` .
Upon return, it contains the value of the Hessian elements
that are possibly non-zero (and will have the same size as *row_out* ).
{xrst_toc_hidden
   example/private/laplace_obj_hes.cpp
}
Example
*******
The file :ref:`laplace_obj_hes.cpp-name` contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

{xrst_end laplace_obj_hes}
*/


// ----------------------------------------------------------------------------
// laplace_obj_hes
void cppad_mixed::laplace_obj_hes(
   const d_vector&          fixed_vec   ,
   const d_vector&          random_vec  ,
   const d_vector&          weight      ,
   CppAD::vector<size_t>&   row_out     ,
   CppAD::vector<size_t>&   col_out     ,
   d_vector&                val_out     )
{  assert( init_laplace_obj_hes_done_ );
   //
   assert( n_fixed_  == fixed_vec.size() );
   assert( n_random_ == random_vec.size() );
   assert( A_rcv_.nr() + 1 == weight.size() );
   //
   size_t nnz = laplace_obj_hes_.subset.nnz();
   if( nnz == 0 )
   {  // sparse Hessian has no entries
      assert( row_out.size() == 0 );
      assert( col_out.size() == 0 );
      assert( val_out.size() == 0 );
      return;
   }
   // check if this is first call
   if( row_out.size() == 0 )
   {  assert( col_out.size() == 0 );
      row_out = laplace_obj_hes_.subset.row();
      col_out = laplace_obj_hes_.subset.col();
      val_out.resize(nnz);
   }
# ifndef NDEBUG
   else
   {  for(size_t k = 0; k < nnz; k++)
      {  assert( row_out[k] == laplace_obj_hes_.subset.row()[k] );
         assert( col_out[k] == laplace_obj_hes_.subset.col()[k] );
      }
   }
# endif
   assert( row_out.size() == nnz );
   assert( col_out.size() == nnz );
   assert( val_out.size() == nnz );
   //
   // beta
   const d_vector& beta(fixed_vec);
   //
   // theta_u
   d_vector theta_u(n_fixed_ + n_random_);
   pack(fixed_vec, random_vec, theta_u);
   //
   // set dynamic parameters in laplace_obj_fun_
   laplace_obj_fun_.new_dynamic(theta_u);
   //
   // val_out
   sparse_rc   not_used_pattern;
   std::string not_used_coloring;
   laplace_obj_fun_.sparse_hes(
      beta                   ,
      weight                 ,
      laplace_obj_hes_.subset ,
      not_used_pattern       ,
      not_used_coloring      ,
      laplace_obj_hes_.work
   );
   for(size_t k = 0; k < nnz; k++)
      val_out[k] = laplace_obj_hes_.subset.val()[k];
   if( CppAD::hasnan( val_out ) ) CppAD::mixed::exception(
      "laplace_obj_hes", "result has a nan"
   );
   //
   return;
}
