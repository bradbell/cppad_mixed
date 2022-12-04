// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>

/*
{xrst_begin fix_like_hes}

Hessian of Fixed Likelihood
###########################

Syntax
******

| *mixed_object* . ``fix_like_hes`` (
| |tab| *fixed_vec* , *weight* , *row_out* , *col_out* , *val_out*
| )

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

fix_likelihood
**************
In the special case where :ref:`fix_likelihood-name` returns the empty vector,
vectors *row_out* , *col_out* , and *val_out* are empty.

fixed_vec
*********
This argument has prototype

   ``const CppAD::vector<double>&`` *fixed_vec*

It specifies the value of the
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>`
vector :math:`\theta` at which the Hessian is evaluated.

weight
******
This argument has prototype

   ``const CppAD::vector<double>&`` *weight*

It specifies the value of the weights for the
components of the
:ref:`problem@Negative Log-Density Vector`
corresponding to the :ref:`fix_likelihood-name` .
It has the same size as the corresponding return value
:ref:`fix_likelihood@vec` .

Hessian
*******
We use :math:`w` to denote the vector corresponding to *weight*
and :math:`v( \theta )` to denote the function corresponding th
the negative log-density vector.
The Hessian is for the function

.. math::

   \sum_{i} w_i v_i ( \theta )

row_out
*******
This argument has prototype

   ``CppAD::vector<size_t>&`` *row_out*

If the input size of this array is non-zero,
the entire vector must be the same
as for a previous call to ``fix_like_hes`` .
If it's input size is zero,
upon return it contains the row indices for the Hessian elements
that are possibly non-zero.

col_out
*******
This argument has prototype

   ``CppAD::vector<size_t>&`` *col_out*

If the input size of this array is non-zero,
the entire vector must be the same as for
a previous call to ``fix_like_hes`` .
If it's input size is zero,
upon return it contains the column indices for the Hessian elements
that are possibly non-zero (and will have the same size as *row_out* ).

val_out
*******
This argument has prototype

   ``CppAD::vector<double>&`` *val_out*

If the input size of this array is non-zero, it must have the same size
as for a previous call to ``fix_like_hes`` .
Upon return, it contains the value of the Hessian elements
that are possibly non-zero (and will have the same size as *row_out* ).
{xrst_toc_hidden
   example/private/fix_like_hes.cpp
}
Example
*******
The file :ref:`fix_like_hes.cpp-name` contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

{xrst_end fix_like_hes}
*/


void cppad_mixed::fix_like_hes(
   const d_vector&        fixed_vec   ,
   const d_vector&        weight      ,
   CppAD::vector<size_t>& row_out     ,
   CppAD::vector<size_t>& col_out     ,
   d_vector&              val_out     )
{  assert( init_fix_like_done_ );
   //
   size_t nnz = fix_like_hes_.subset.nnz();
   if( nnz == 0 )
   {  // sparse Hessian has no entries
      assert( row_out.size() == 0 );
      assert( col_out.size() == 0 );
      assert( val_out.size() == 0 );
      return;
   }
   if( row_out.size() == 0 )
   {  assert( col_out.size() == 0 );
      row_out = fix_like_hes_.subset.row();
      col_out = fix_like_hes_.subset.col();
      val_out.resize(nnz);
   }
# ifndef NDEBUG
   else
   {  for(size_t k = 0; k < nnz; k++)
      {  assert( row_out[k] == fix_like_hes_.subset.row()[k] );
         assert( col_out[k] == fix_like_hes_.subset.col()[k] );
      }
   }
# endif
   assert( row_out.size() == nnz );
   assert( col_out.size() == nnz );
   assert( val_out.size() == nnz );
   //
   sparse_rc   not_used_pattern;
   std::string not_used_coloring;
   fix_like_fun_.sparse_hes(
      fixed_vec            ,
      weight               ,
      fix_like_hes_.subset ,
      not_used_pattern     ,
      not_used_coloring    ,
      fix_like_hes_.work
   );
   for(size_t k = 0; k < nnz; k++)
      val_out[k] = fix_like_hes_.subset.val()[k];
   //
   //
   if( CppAD::hasnan( val_out ) ) throw CppAD::mixed::exception(
      "fix_like_hes", "result has a nan"
   );
   return;
}
