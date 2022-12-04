// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/ipopt_random.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
{xrst_begin ipopt_random_get_bounds_info}

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
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

Prototype
*********
{xrst_spell_off}
{xrst_code cpp} */
bool ipopt_random::get_bounds_info(
      Index       n        ,   // in
      Number*     x_l      ,   // out
      Number*     x_u      ,   // out
      Index       m        ,   // in
      Number*     g_l      ,   // out
      Number*     g_u      )   // out
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_random_get_bounds_info}
*/
{
   assert( size_t(n) == n_random_ );
   assert( m == 0 );
   //
   for(size_t j = 0; j < n_random_; j++)
   {  // map infinity to crazy value required by ipopt
      if( random_lower_[j] == - std::numeric_limits<double>::infinity() )
         x_l[j] = nlp_lower_bound_inf_;
      else
         x_l[j] = random_lower_[j];
      //
      if( random_upper_[j] == std::numeric_limits<double>::infinity() )
         x_u[j] = nlp_upper_bound_inf_;
      else
         x_u[j] = random_upper_[j];
   }
   //
   return true;
}
} } // END_CPPAD_MIXED_NAMESPACE
