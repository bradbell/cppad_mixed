// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin ldlt_cholmod_logdet}
{xrst_spell
   determinant
   logdet
}

Compute Log Determinant for Current Factor
##########################################

Syntax
******

   *logdet* = *ldlt_obj* . ``logdet`` ( *negative* )

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

ldlt_obj
********
This object has prototype

   ``const CppAD::mixed::ldlt_cholmod`` *ldlt_obj*

In addition, it must have a previous call to
:ref:`ldlt_cholmod_update-name` .

negative
********
The input value of this argument does no matter.
Upon return, it is the number of elements of
:ref:`ldlt_cholmod@Factorization@D`
that are less than zero.

logdet
******
This return value
is the log of the absolute value of the determinant corresponding
to the previous call to :ref:`ldlt_cholmod_update-name` .
If the matrix is singular, *logdet* is
minus infinity.

Example
*******
The file :ref:`ldlt_cholmod.cpp<ldlt_cholmod.cpp@logdet>` contains an
example and test that uses this function.

{xrst_end ldlt_cholmod_logdet}
*/

# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cassert>
# include <cmath>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN_PROTOTYPE
double ldlt_cholmod::logdet(size_t& negative) const
// END_PROTOTYPE
{  assert( update_called_ );
   //
   // factorization P H P' = L D L'
   int*    L_p  = (int *) factor_->p;
   double* L_x  = (double *) factor_->x;
# ifndef NDEBUG
   int*    L_i  = (int *) factor_->i;
# endif
   negative        = 0;
   bool has_zero   = false;
   double logdet_H = 0.0;
   for(size_t j = 0; j < nrow_; j++)
   {  // first element for each column is always the diagonal element
      assert( size_t( L_i [ L_p[j] ] ) == j );
      // j-th element on diagonal of D in factorization
      double dj = L_x[ L_p[j] ];
      has_zero |= dj == 0.0;
      if( dj < 0.0 )
         negative++;
      if( ! has_zero )
         logdet_H += std::log( std::fabs(dj) );
   }
   if( has_zero )
      return - std::numeric_limits<double>::infinity();
   //
   return logdet_H;
}

} } // END_CPPAD_MIXED_NAMESPACE
