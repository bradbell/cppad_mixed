// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/ipopt_random.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
{xrst_begin ipopt_random_eval_grad_f}

Compute Gradient of the Objective
#################################

Syntax
******
*ok* = ``eval_grad_f`` ( *n* , *x* , *new_x* , *grad_f* )

n
*
is the number of variables in the problem (dimension of x).

x
*
is the value for the primal variables at which the gradient
:math:`\nabla f(x)` is computed (has size *n* ).

new_x
*****
if true, no Ipopt evaluation method was previous called with the same
value for *x* .

grad_f
******
is set to the value for the gradient :math:`\nabla f(x)`
(has size *m* ).

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
bool ipopt_random::eval_grad_f(
   Index           n         ,  // in
   const Number*   x         ,  // in
   bool            new_x     ,  // in
   Number*         grad_f    )  // out
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_random_eval_grad_f}
*/
{  try
   {  try_eval_grad_f(n, x, new_x, grad_f);
   }
   catch(const std::exception& e)
   {  error_message_ = "ipopt_random::eval_grad_f: std::exception: ";
      for(size_t j = 0; j < n_random_; j++)
         error_random_[j] = x[j];
      return false;
   }
   catch(const CppAD::mixed::exception& e)
   {  error_message_ = e.message("ipopt_random::eval_grad_f");
      for(size_t j = 0; j < n_random_; j++)
         error_random_[j] = x[j];
      return false;
   }
   return true;
}
void ipopt_random::try_eval_grad_f(
   Index           n         ,  // in
   const Number*   x         ,  // in
   bool            new_x     ,  // in
   Number*         grad_f    )  // out
{  assert( size_t(n) == n_random_ );
   assert( mixed_object_.ran_like_fun_.Domain() == n_fixed_ + n_random_ );
   assert( mixed_object_.ran_like_fun_.Range()  == 1 );
   //
   if( new_x )
   {  // set the zero order Taylor coefficients in
      // mixed_object_.ran_like_fun_
      Number obj_value;
      eval_f(n, x, new_x, obj_value);
   }
   // compute the gradient w.r.t fixed and random effects
   d_vector w(1);
   w[0] = 1.0;
   d_vector dw = mixed_object_.ran_like_fun_.Reverse(1, w);
   if( CppAD::hasnan( dw ) ) throw CppAD::mixed::exception(
      "", "gradient has nan"
   );
   //
   // return gradient w.r.t random effects
   for(size_t j = 0; j < n_random_; j++)
      grad_f[j] = dw[ n_fixed_ + j ];
   return;
}
} } // END_CPPAD_MIXED_NAMESPACE
