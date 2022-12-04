// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/ipopt_random.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
{xrst_begin ipopt_random_eval_g}

Compute Value of Constraint Functions
#####################################

Syntax
******
*ok* = ``eval_g`` ( *n* , *x* , *new_x* , *m* , *g* )

n
*
is the number of variables in the problem (dimension of x).

x
*
is the value for the primal variables at which the constraints
:math:`g(x)` is computed (has size *n* ).

new_x
*****
if true, no Ipopt evaluation method was previous called with the same
value for *x* .

m
*
is the number of constraints in the problem (dimension of g(x)).

g
*
is set to the value for the constraint functions (has size *m* ).

ok
**
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

mixed_object\_.ran_like_fun\_
*****************************
if *new_x* is true,
after this call, the zero order Taylor coefficients in this function
will corresponding to the value of *fixed_vec_* and
the random effects in *x* .

Prototype
*********
{xrst_spell_off}
{xrst_code cpp} */
bool ipopt_random::eval_g(
   Index           n        ,  // in
   const Number*   x        ,  // in
   bool            new_x    ,  // in
   Index           m        ,  // in
   Number*         g        )  // out
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_random_eval_g}
*/
{  assert( size_t(n) == n_random_ );
   assert( m == 0 );
   //
   if( new_x )
   {  // set the zero order Taylor coefficients in
      // mixed_object_.ran_like_fun_
      Number obj_value;
      eval_f(n, x, new_x, obj_value);
   }
   //
   return true;
}
} } // END_CPPAD_MIXED_NAMESPACE
