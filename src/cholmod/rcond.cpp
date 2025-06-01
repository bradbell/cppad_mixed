// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2025 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin ldlt_cholmod_rcond dev}

Reciprocal of Condition Number for D
####################################

Syntax
******
*rcond* = *ldlt_obj* . ``rcond`` ( )

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE
   // END_PROTOTYPE
}

Private
*******
The ``ldlt_cholmod`` class is an
:ref:`implementation detail<ldlt_cholmod@Private>` and not part of the
CppAD Mixed user API.

ldlt_obj
********
The object *ldlt_obj*
must have a previous call to :ref:`ldlt_cholmod_update-name` .

rcond
*****
This return value *rcond*
is the reciprocal of the condition number for the diagonal matrix *D*
in the factorization.
In other words, it is the minimum absolute entry in *D* divided
by the maximum absolute entry in *D* .
If the matrix *D* is singular, or any entry in *D* nan or infinite,
*rcond* is zero.

Example
*******
The file :ref:`ldlt_cholmod.cpp<ldlt_cholmod.cpp@rcond>` contains an
example and test that uses this function.

{xrst_end ldlt_cholmod_rcond}
*/

# include <cppad/mixed/ldlt_cholmod.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN_PROTOTYPE
double ldlt_cholmod::rcond(void) const
// END_PROTOTYPE
{  assert( update_called_ );
   //
   // factorization P H P' = L D L'
   int*    L_p  = (int *) factor_->p;
   double* L_x  = (double *) factor_->x;
# ifndef NDEBUG
   int*    L_i  = (int *) factor_->i;
# endif
   double max_abs  = 0.0;
   double min_abs  = CppAD::numeric_limits<double>::infinity();
   for(size_t j = 0; j < nrow_; j++)
   {  // first element for each column is always the diagonal element
      assert( size_t( L_i [ L_p[j] ] ) == j );
      // j-th element on diagonal of D in factorization
      double dj  = L_x[ L_p[j] ];
      double abs = fabs( dj );
      if( isnan( abs ) )
         abs = 0.0;
      max_abs    = std::max( abs, max_abs);
      min_abs    = std::min( abs, min_abs);
   }
   //
   // rcond
   if( min_abs == 0.0 )
      return 0.0;
   if( min_abs == CppAD::numeric_limits<double>::infinity() )
      return 0.0;
   if( max_abs == CppAD::numeric_limits<double>::infinity() )
      return 0.0;
   return min_abs / max_abs;
}
} } // END_CPPAD_MIXED_NAMESPACE
