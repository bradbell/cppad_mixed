# ifndef CPPAD_MIXED_ONE_DIM_DERIVATIVE_CHK_HPP
# define CPPAD_MIXED_ONE_DIM_DERIVATIVE_CHK_HPP

// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin one_dim_derivative_chk}
{xrst_spell
   chk
   dfdx
   rel
   tol
}

One Dimensional Finite Difference Derivative Check
##################################################

Syntax
******

| *result* = ``one_dim_derivative_chk`` (
| |tab| *obj* , *x_lower* , *x_upper* , *x_mid* , *f_mid* , *dfdx* , *rel_tol*
| )

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE
   // END_PROTOTYPE
}

obj
***
This is a object such that the syntax

   *ok* = *obj* . ``one_dim_function`` ( *x_in* , *f_out* )

evaluates the function.
The argument *x_in* has prototype

   ``double x``

and is the point at which we are evaluating the function.
The argument *f_out* has prototype

   ``d_vector& f_out``

Its input value does not matter and upon return it is the
value of the function at *x_in* .
The return value *ok* has prototype

   ``bool`` *ok*

If it is true, the computation of *f_out* completed.
Otherwise, an error was detected and the computation was aborted.

x_lower
*******
This is a lower limit for the argument *x_in* to *fun* .
If it is less than infinity, it is used to get an approximation
for the scale of the argument to the function.

x_upper
*******
This is a lower limit for the argument *x_in* to *fun* .
If it is greater than - infinity, it is used to get an approximation
for the scale of the argument to the function.

x_mid
*****
Is the argument value at which the derivative is checked.
It is the mid-point of a central difference approximation for the derivative
(unless the bounds *x_lower* , *x_upper* prevent this.)

f_mid
*****
Is the value of the function at *x_mid* .

dfdx
****
Is the value of the derivative that we are checking with a central difference
approximation.

rel_tol
*******
Is an acceptable relative tolerance for the difference between
*dfdx* and its finite difference approximation.

result
******
The structure *result* [ *i* ] contains
the following information about the *i*-th component of the function:

rel_err
=======
The smallest relative error that ``one_dim_derivative_check`` found,
for the *i*-th component of the function.
If it is less than or equal *rel_tol* the derivative check passed.
Searching for a step size with a smaller relative error stops as soon
as the derivative check passes for all components of the function.
If *result* [ *i* ]. ``rel_err`` is infinity,
then we were not able to evaluate the function.

step
====
The step size corresponding to *result* [ *i* ]. ``rel_err`` .

f_minus
=======
Is the value of f(x) at

   *x_minus* = ``max`` ( *x_mid* ``-`` *step* , *x_minus* )

f_plus
======
Is the value of f(x) at

   *x_plus* = ``min`` ( *x_mid* + *step* , *x_upper* )

Derivative Approximation
========================
The corresponding derivative approximation is

   ( *f_plus* ``-`` *f_minus* ) / ( *x_plus* ``-`` *x_minus* )

{xrst_end one_dim_derivative_chk}
*/
# include <cppad/mixed/typedef.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
// BEGIN_PROTOTYPE
struct one_dim_derivative_result {
   double rel_err;
   double step;
   double f_minus;
   double f_plus;
};
template <class Object>
CppAD::vector<one_dim_derivative_result> one_dim_derivative_chk(
   Object&    obj          ,
   double     x_lower      ,
   double     x_upper      ,
   double     x            ,
   d_vector   f            ,
   d_vector   dfdx         ,
   double     rel_tol      )
// END_PROTOTYPE
{  // infinity
   double infinity = std::numeric_limits<double>::infinity();
   //
   // nan
   double nan      = std::numeric_limits<double>::quiet_NaN();
   //
   // sqrt_eps
   double sqrt_eps = std::sqrt( std::numeric_limits<double>::epsilon() );
   //
   // m
   size_t m = f.size();
   assert( m == dfdx.size() );
   //
   // x_max_abs
   double x_max_abs = std::fabs(x);
   if( std::fabs(x_upper) < infinity )
      x_max_abs = std::max(x_max_abs, std::fabs(x_upper));
   if( std::fabs(x_lower) < infinity )
      x_max_abs = std::max(x_max_abs, std::fabs(x_lower));
   //
   // log_max_rel_step
   // log of the maximum relative step to try
   double log_max_rel_step = std::log(0.1);
   if( x_max_abs == 0.0 )
      log_max_rel_step = std::log(1e+5);
   //
   // log_min_rel_step
   // log of the minimum relative step to try
   double log_min_rel_step = std::log(1e-10);
   //
   // n_try
   // number of relative steps to try
   size_t n_try = 5;
   if( x_max_abs == 0.0 )
      n_try = 7;
   //
   // log_diff
   // difference of log of relative step between trys
   double log_diff = (log_max_rel_step - log_min_rel_step)/double(n_try - 1);
   //
   // result
   // initialize
   CppAD::vector<one_dim_derivative_result> result(m);
   for(size_t i = 0; i < m; ++i)
   {  result[i].rel_err   = infinity;
      result[i].step      = nan;
      result[i].f_minus   = nan;
      result[i].f_plus    = nan;
   }
   //
   // loop over finite difference step sizes
   size_t i_try       = 0;        // index for this step size
   double rel_err_max = infinity; // maximum in result
   d_vector f_plus(m), f_minus(m);
   while( i_try < n_try && rel_err_max > rel_tol )
   {  // rel_step
      double log_try  = log_min_rel_step + log_diff * double(i_try);
      double rel_step = std::exp(log_try);
      //
      // step
      double step     = rel_step;
      if( x_max_abs != 0.0 )
         step = rel_step * x_max_abs;
      //
      // x_plus
      double x_plus = std::min(x + step, x_upper);
      //
      // x_minus
      double x_minus = std::max(x - step, x_lower);
      //
      // f_plus
      bool ok = obj.one_dim_function(x_plus, f_plus);
      //
      // f_minus
      ok = ok && obj.one_dim_function(x_minus, f_minus);
      //
      // rel_err_max
      rel_err_max = 0.0;
      for(size_t i = 0; i < m; ++i)
      {  // f_max_abs
         double f_max_abs = std::fabs(f[i]);
         f_max_abs        = std::max(f_max_abs, std::fabs(f_plus[i]));
         f_max_abs        = std::max(f_max_abs, std::fabs(f_minus[i]));
         //
         // apx_dfdx
         // finite difference approximation for derivative
         double apx_dfdx = (f_plus[i] - f_minus[i]) / (x_plus - x_minus);
         //
         // rel_err
         double rel_err = std::fabs(dfdx[i] - apx_dfdx);
         double den     = std::fabs(dfdx[i]) + std::fabs(apx_dfdx);
         den           += sqrt_eps * f_max_abs;
         if( den > 0.0 )
            rel_err = rel_err / den;
         //
         // best so far ?
         if( ok && 1.1 * rel_err < result[i].rel_err )
         {  result[i].rel_err  = rel_err;
            result[i].step     = step;
            result[i].f_minus  = f_minus[i];
            result[i].f_plus   = f_plus[i];
         }
         // rel_err_max
         rel_err_max = std::max(rel_err_max, result[i].rel_err);
      }
      //
      // next try
      ++i_try;
   }
   return result;
}
} } // END_CPPAD_MIXED_NAMESPACE

# endif
