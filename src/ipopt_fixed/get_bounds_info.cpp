// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/ipopt_fixed.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
{xrst_begin ipopt_fixed_get_bounds_info}

Return Optimization Bounds
##########################

Syntax
******
*ok* = ``get_bounds_info`` ( *n* , *x_l* , *x_u* , *m* , *g_l* , *g_u* )

n
*
is the number of variables in the problem (dimension of x).

x_l
***
set to the lower bounds for *x* (has size *n* ).

x_u
***
set to the upper bounds for *x* (has size *n* ).

m
*
is the number of constraints in the problem (dimension of g(x)).

g_l
***
set to the lower bounds for *g* ( *x* ) (has size *m* ).

g_u
***
set to the upper bounds for *g* ( *x* ) (has size *m* ).

ok
**
If set to false, the optimization will terminate with status set to
:ref:`ipopt_fixed_finalize_solution@status@USER_REQUESTED_STOP` .

Prototype
*********
{xrst_spell_off}
{xrst_code cpp} */
bool ipopt_fixed::get_bounds_info(
      Index       n        ,   // in
      Number*     x_l      ,   // out
      Number*     x_u      ,   // out
      Index       m        ,   // in
      Number*     g_l      ,   // out
      Number*     g_u      )   // out
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_fixed_get_bounds_info}
*/
{  assert( adaptive_called_ );
   double inf = std::numeric_limits<double>::infinity();
   //
   assert( n > 0 );
   assert( size_t(n) == n_fixed_ + fix_likelihood_nabs_ );
   assert( m >= 0 );
   assert( size_t(m) == 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_ );

   // fixed effects
   for(size_t j = 0; j < n_fixed_; j++)
   {  // map infinity to crazy value required by ipopt
      if( fixed_lower_[j] == - inf )
         x_l[j] = nlp_lower_bound_inf_;
      else
         x_l[j] = scale_x_[j] * fixed_lower_[j];
      //
      if( fixed_upper_[j] == inf )
         x_u[j] = nlp_upper_bound_inf_;
      else
         x_u[j] = scale_x_[j] * fixed_upper_[j];
      if( fixed_lower_[j] == fixed_upper_[j] )
         x_u[j] = x_l[j];
   }
   // auxillary varibles for absolute value terms
   for(size_t j = 0; j < fix_likelihood_nabs_; j++)
   {  x_l[n_fixed_ + j] = nlp_lower_bound_inf_;
      x_u[n_fixed_ + j] = nlp_upper_bound_inf_;
   }
   //
   // constraints for absolute value terms
   for(size_t j = 0; j < 2 * fix_likelihood_nabs_; j++)
   {  g_l[j] = 0.0;
      g_u[j] = nlp_upper_bound_inf_;
   }
   //
   // fixed constraints
   for(size_t j = 0; j < n_fix_con_; j++)
   {  size_t i = 2 * fix_likelihood_nabs_ + j;
      // g_l
      if( fix_constraint_lower_[j] == - inf )
         g_l[i] = nlp_lower_bound_inf_;
      else
         g_l[i] = scale_g_[i] * fix_constraint_lower_[j];
      // g_u
      if( fix_constraint_upper_[j] == inf )
         g_u[i] = nlp_upper_bound_inf_;
      else
         g_u[i] = scale_g_[i] * fix_constraint_upper_[j];
      if( fix_constraint_lower_[j] == fix_constraint_upper_[j] )
         g_u[i] = g_l[i];
   }
   //
   // random constraints
   for(size_t j = 0; j < n_ran_con_; j++)
   {  g_l[2 * fix_likelihood_nabs_ + n_fix_con_ + j] = 0.0;
      g_u[2 * fix_likelihood_nabs_ + n_fix_con_ + j] = 0.0;
   }
   //
   return true;
}
} } // END_CPPAD_MIXED_NAMESPACE
