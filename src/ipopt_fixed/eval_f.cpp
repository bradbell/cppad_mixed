// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/ipopt_fixed.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
{xrst_begin ipopt_fixed_eval_f dev}

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
If set to false, the optimization will terminate with status set to
:ref:`ipopt_fixed_finalize_solution@status@USER_REQUESTED_STOP` .

Prototype
*********
{xrst_spell_off}
{xrst_code cpp} */
bool ipopt_fixed::eval_f(
   Index           n         ,  // in
   const Number*   x         ,  // in
   bool            new_x     ,  // in
   Number&         obj_value )  // out
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_fixed_eval_f}
*/
{  for(Index j = 0; j < n; ++j)
      x_tmp_[j] = scale_x_[j] * x[j];
   if( abort_on_eval_error_ )
   {  try_eval_f(n, x_tmp_, new_x, obj_value);
   }
   else
   {  try
      {  try_eval_f(n, x_tmp_, new_x, obj_value);
      }
      catch(const std::exception& e)
      {  error_message_ = "ipopt_fixed::eval_f: std::exception: ";
         for(size_t j = 0; j < n_fixed_; j++)
            error_fixed_[j] = x[j];
         return false;
      }
      catch(const CppAD::mixed::exception& e)
      {  error_message_ = e.message("ipopt_fixed::eval_f");
         for(size_t j = 0; j < n_fixed_; j++)
            error_fixed_[j] = x[j];
         return false;
      }
   }
   obj_value *= scale_f_;
   return true;
}
void ipopt_fixed::try_eval_f(
   Index           n         ,  // in
   const d_vector& x         ,  // in
   bool            new_x     ,  // in
   Number&         obj_value )  // out
{
   assert( n > 0 && size_t(n) == n_fixed_ + fix_likelihood_nabs_ );
   //
   // value of fixed effects corresponding to this x
   for(size_t j = 0; j < n_fixed_; j++)
      fixed_tmp_[j] = double( x[j] );
   //
   // random part of objective
   double H = Number( 0.0 );
   if( n_random_ > 0 )
   {  //
      // compute the optimal random effects corresponding to fixed effects
      if( new_x )
         new_random(fixed_tmp_);
      H = mixed_object_.ran_obj_eval(fixed_tmp_, random_cur_);
   }
   obj_value = Number(H);
   if( fix_likelihood_vec_tmp_.size() == 0 )
      assert( mixed_object_.fix_like_eval(fixed_tmp_).size() == 0 );
   else
   {
      // fixed part of objective
      // (2DO: cache fix_likelihood_vec_tmp_ for eval_g with same x)
      fix_likelihood_vec_tmp_ = mixed_object_.fix_like_eval(fixed_tmp_);
      //
      // only include smooth part of prior in objective
      obj_value += Number( fix_likelihood_vec_tmp_[0] );
      //
      // auxillary variable with index j is constrainted to be
      // greater than absolute value of 1+j component of fixed likelihood
      for(size_t j = 0; j < fix_likelihood_nabs_; j++)
         obj_value += x[n_fixed_ + j];
   }
   if( CppAD::isnan(obj_value) ) throw CppAD::mixed::exception(
      "", "objective is nan"
   );
   return;
}
} } // END_CPPAD_MIXED_NAMESPACE
