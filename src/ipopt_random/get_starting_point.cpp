// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/ipopt_random.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
{xrst_begin ipopt_random_get_starting_point dev}

Return Initial Values Where Optimization is Started
###################################################

Syntax
******

| *ok* = ``get_starting_point`` (
| |tab| *n* , *init_x* , *x* , *init_z* , *z_L* , *z_U* , *m* , *init_lambda* , *lambda*
| )

n
*
is the number of variables in the problem (dimension of x).

init_x
******
assumed true which means the ipopt options specify that the this routine
will provide an initial value for *x* .

x
*
if *init_x* is true,
set to the initial value for the primal variables (has size *n* ).

init_z
******
assumes *init_z* is false.
If it were true, the ipopt options specify that the this routine
will provide an initial value for *x* upper and lower bound
multipliers.

z_L
***
if *init_z* is true,
set to the initial value for the lower bound multipliers (has size *n* ).

z_U
***
if *init_z* is true,
set to the initial value for the upper bound multipliers (has size *n* ).

init_lambda
***********
assumes *init_lambda* is false.
If it were true, the ipopt options specify that the this routine
will provide an initial value for *g* ( *x* ) upper and lower bound
multipliers.

lambda
******
if *init_lambda* is true,
set to the initial value for the *g* ( *x* ) multipliers
(has size *m* ).

ok
**
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

Prototype
*********
{xrst_spell_off}
{xrst_code cpp} */
bool ipopt_random::get_starting_point(
   Index           n            ,  // in
   bool            init_x       ,  // in
   Number*         x            ,  // out
   bool            init_z       ,  // in
   Number*         z_L          ,  // out
   Number*         z_U          ,  // out
   Index           m            ,  // in
   bool            init_lambda  ,  // in
   Number*         lambda       )  // out
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_random_get_starting_point}
*/
{
   assert( init_x == true );
   assert( init_z == false );
   assert( init_lambda == false );
   assert( size_t(n) == n_random_ );
   assert( m == 0 );
   //
   // initial value for random effects
   for(size_t j = 0; j < n_random_; j++)
      x[j] = random_in_[j];

   return true;
}
} } // END_CPPAD_MIXED_NAMESPACE
