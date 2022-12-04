// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/ipopt_random.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
{xrst_begin ipopt_random_eval_f}

Compute Value of Objective
##########################

Syntax
******
*ok* = ``eval_f`` ( *n* , *x* , *new_x* , *obj_value* )

n
*
is the number of variables in the problem (dimension of x).

x
*
is the value for the primal variables at which the objective
f(x) is computed (has size *n* ).

new_x
*****
if true, no Ipopt evaluation method was previous called with the same
value for *x* .

obj_val
*******
set to the initial value of the objective function f(x).

ok
**
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

objective_current\_
*******************
After this call, ``objective_current_`` will be the
value of the objective corresponding to *x* .
Note that if *new_x* is false,
this value does not change.

mixed_object\_.ran_like_fun\_
*****************************
After this call, the zero order Taylor coefficients in this function
will corresponding to the value of *fixed_vec_* and
the random effects in *x* .

Prototype
*********
{xrst_spell_off}
{xrst_code cpp} */
bool ipopt_random::eval_f(
   Index           n         ,  // in
   const Number*   x         ,  // in
   bool            new_x     ,  // in
   Number&         obj_value )  // out
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_random_eval_f}
*/
{  try
   {  try_eval_f(n, x, new_x, obj_value);
   }
   catch(const std::exception& e)
   {  error_message_ = "ipopt_random::eval_f: std::exception: ";
      for(size_t j = 0; j < n_random_; j++)
         error_random_[j] = x[j];
      return false;
   }
   catch(const CppAD::mixed::exception& e)
   {  error_message_ = e.message("ipopt_random::eval_f");
      for(size_t j = 0; j < n_random_; j++)
         error_random_[j] = x[j];
      return false;
   }
   return true;
}
void ipopt_random::try_eval_f(
   Index           n         ,  // in
   const Number*   x         ,  // in
   bool            new_x     ,  // in
   Number&         obj_value )  // out
{  assert( size_t(n) == n_random_ );
   assert( mixed_object_.ran_like_fun_.Domain() == n_fixed_ + n_random_ );
   assert( mixed_object_.ran_like_fun_.Range()  == 1 );
   //
   if( new_x )
   {  // random effects as a vector
      d_vector random_vec(n_random_);
      for(size_t j = 0; j < n_random_; j++)
         random_vec[j] = x[j];
      //
      // pack both the fixed and random effects into one vector
      d_vector both_vec(n_fixed_ + n_random_);
      mixed_object_.pack(fixed_vec_, random_vec, both_vec);
      //
      // compute the log-density vector
      d_vector vec = mixed_object_.ran_like_fun_.Forward(0, both_vec);
      assert( vec.size() == 1 );
      if( CppAD::hasnan( vec ) ) throw CppAD::mixed::exception(
         "", "objective has a nan"
      );
      // store for re-use
      objective_current_ = vec[0];
   }
   //
   obj_value = objective_current_;
   return;
}
} } // END_CPPAD_MIXED_NAMESPACE
