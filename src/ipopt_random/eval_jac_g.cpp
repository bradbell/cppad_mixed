// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/ipopt_random.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
{xrst_begin ipopt_random_eval_jac_g}
{xrst_spell
   nele
}

Compute Jacobian of Constraint Functions
########################################

Syntax
******

| *ok* = ``eval_jac_g`` (
| |tab| *n* , *x* , *new_x* , *m* , *nele_jac* , *iRow* , *jCol* , *values*
| )

n
*
is the number of variables in the problem (dimension of x).

x
*
is the value for the primal variables at which the Jacobian
of the constraints :math:`\nabla g(x)` is computed (has size *n* ).

new_x
*****
if true, no Ipopt evaluation method was previous called with the same
value for *x* .

m
*
is the number of constraints in the problem (dimension of g(x)).

nele_jac
********
is the number of non-zero elements in the Jacobian of *g* ( *x* ) ; i.e.,
the same as
:ref:`ipopt_random_get_nlp_info@nnz_jac_g` .

iRow
****
If *values* is ``NULL`` ,
*iRow* has size *nele_jac* and is set to the
row indices for the non-zero entries in the Jacobian of the constraints
:math:`g_x (x)`.

jCol
****
If *values* is ``NULL`` ,
*jCol* has size *nele_jac* and is set to the
column indices for the non-zero entries in the Jacobian of the constraints
:math:`g_x (x)`.

values
******
If *values* is not ``NULL`` ,
it has size *nele_jac* and *values* [ *k* ]
is set to the value of element of the Jacobian :math:`g_x (x)`
with row index *iRow* [ *k* ]
and column index *jCol* [ *k* ] .

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
bool ipopt_random::eval_jac_g(
   Index           n        ,  // in
   const Number*   x        ,  // in
   bool            new_x    ,  // in
   Index           m        ,  // in
   Index           nele_jac ,  // in
   Index*          iRow     ,  // out
   Index*          jCol     ,  // out
   Number*         values   )  // out
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_random_eval_jac_g}
*/
{  try
   {  try_eval_jac_g(n, x, new_x, m, nele_jac, iRow, jCol, values);
   }
   catch(const std::exception& e)
   {  error_message_ = "ipopt_random::eval_jac_g: std::exception: ";
      for(size_t j = 0; j < n_random_; j++)
         error_random_[j] = x[j];
      return false;
   }
   catch(const CppAD::mixed::exception& e)
   {  error_message_ = e.message("ipopt_random::eval_jac_g");
      for(size_t j = 0; j < n_random_; j++)
         error_random_[j] = x[j];
      return false;
   }
   return true;
}
void ipopt_random::try_eval_jac_g(
   Index           n        ,  // in
   const Number*   x        ,  // in
   bool            new_x    ,  // in
   Index           m        ,  // in
   Index           nele_jac ,  // in
   Index*          iRow     ,  // out
   Index*          jCol     ,  // out
   Number*         values   )  // out
{  assert( size_t(n) == n_random_ );
   assert( m == 0 );
   assert( nele_jac == 0 );
   if( values == NULL )
   {  assert( ! new_x );
      return;
   }
   //
   if( new_x )
   {  // set the zero order Taylor coefficients in
      // mixed_object_.ran_like_fun_
      Number obj_value;
      eval_f(n, x, new_x, obj_value);
   }
   return;
}
} } // END_CPPAD_MIXED_NAMESPACE
