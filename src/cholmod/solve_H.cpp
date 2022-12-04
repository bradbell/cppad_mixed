// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin ldlt_cholmod_solve_H}

Solve Linear Equations Using Factor
###################################

Syntax
******
*ldlt_obj* . ``solve_H`` ( *row* , *val_in* , *val_out* )

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
This function solves the linear equation
:math:`H x = b` where :math:`H` is the symmetric matrix
that has been factored,
:math:`b` is a known column vector,
and :math:`x` is unknown.

ldlt_obj
********
This object has prototype

   ``const CppAD::mixed::ldlt_cholmod`` *ldlt_obj*

In addition, it must have a previous call to
:ref:`ldlt_cholmod_update-name` .

row
***
This argument
contains all of the rows of column vector :math:`b` that are
non-zero and the rows of the column vector *x*
that are desired.
These values are in strictly increasing order; i.e.,

   *row* [ *k* ] < *row* [ *k* +1]

It follows that *row* . ``size`` () is less than or equal
:ref:`ldlt_cholmod_ctor@nrow_` .

val_in
******
This argument
has the same size as *row* .
It specifies the values in the column vector :math:`b`
for each of the corresponding rows; i.e.,
for *k* = 0 , ..., *row* . ``size`` () ``-1`` ,

   *b* [ *row* [ *k* ] ] = *val_in* [ *k* ]

.

val_out
*******
This argument
has the same size as *row* .
On input, the value of its elements do not matter.
Upon return, it contains the values in the column vector :math:`b`
for each of the corresponding rows; i.e.,
for *k* = 0 , ..., *row* . ``size`` () ``-1`` ,

   *x* [ *row* [ *k* ] ] = *val_out* [ *k* ]

.

Example
*******
The file :ref:`ldlt_cholmod.cpp<ldlt_cholmod.cpp@solve_H>` contains an
example and test that uses this function.

{xrst_end ldlt_cholmod_solve_H}
*/
# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cppad/utility/index_sort.hpp>
# include <cassert>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN_PROTOTYPE
void ldlt_cholmod::solve_H(
   const CppAD::vector<size_t>& row      ,
   const CppAD::vector<double>& val_in   ,
   CppAD::vector<double>&       val_out  )
// END_PROTOTYPE
{  assert( update_called_ );
   assert( row.size() == val_in.size() );
   assert( row.size() == val_out.size() );
   assert( row.size() <= nrow_ );
# ifndef NDEBUG
   for(size_t k = 1; k < row.size(); k++)
      assert( row[k-1] < row[k] && row[k] < nrow_ );
# endif
   //
   assert( rhs_ != CPPAD_NULL  );
   assert( rhs_->nrow == nrow_ );
   assert( rhs_->ncol == 1     );
   double* rhs_x = (double *) rhs_->x;
   //
   assert( rhs_set_ != CPPAD_NULL );
   int* rhs_set_p = (int *) rhs_set_->p;
   int* rhs_set_i = (int *) rhs_set_->i;
   //
   // set non-zero entries in right hand size rhs_
   // (value of rhs_ in other rows does not matter because not in rhs_set)
   for(size_t k = 0; k < row.size(); k++)
      rhs_x[ row[k] ] = val_in[k];
   //
   rhs_set_p[0] = 0;
   rhs_set_p[1] = int( row.size() );
   for(size_t k = 0; k < row.size(); k++)
      rhs_set_i[k] = (int) row[k];

   // solve the linear equation H * sol = rhs
   int sys = CHOLMOD_A;
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
   //
   assert( sol_set_->nrow == nrow_ );
   assert( sol_set_->ncol == 1 );
   assert( sol_set_->xtype == CHOLMOD_PATTERN );
   assert( sol_set_->packed == CHOLMOD_TRUE);
   int* sol_set_p = (int *) sol_set_->p;
   int* sol_set_i = (int *) sol_set_->i;
   //
   assert( sol_ != CPPAD_NULL   );
   assert( sol_->nrow == nrow_ );
   assert( sol_->ncol == 1     );
   double* sol_x = (double *) sol_->x;
   //
   // sort_ is an index sort of sol_set_i
   size_t ni = (size_t) sol_set_p[1];
   key_.resize(ni);
   index_.resize(ni);
   for(size_t ell = 0; ell < ni; ell++)
      key_[ell] = (size_t) sol_set_i[ell];
   CppAD::index_sort(key_, index_);

   // return result values
   size_t k   = 0;
   size_t ell = 0;
   while(ell < ni && k < row.size() )
   {  size_t i = key_[ index_[ell++] ];
      assert( i <= row[k] );
      if( i == row[k] )
      {  val_out[k] = sol_x[i];
         k++;
      }
   }
   assert( k == row.size() );
}

} } // END_CPPAD_MIXED_NAMESPACE
