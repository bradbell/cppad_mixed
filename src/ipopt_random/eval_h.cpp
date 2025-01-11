// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/ipopt_random.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
{xrst_begin ipopt_random_eval_h dev}
{xrst_spell
  hess
  nele
}

Compute the Hessian of the Lagrangian
#####################################

Syntax
******

| *ok* = ``eval_h`` (
| |tab| *n* , *x* , *new_x* , *obj_factor* , *m* , *lambda* , *new_lambda* , *nele_hess* , *iRow* , *jCol* , *values*
| )

Lagrangian
**********
The Lagrangian is defined to be

.. math::

   L(x) = \alpha f(x) + \sum_{i=0}^{m-1} \lambda_i g_i (x)

mixed_object.quasi_fixed\_
**************************
It is assumed that this member variable is false.

n
*
is the number of variables in the problem (dimension of x).

x
*
is the value for the primal variables at which the
Hessian of the Lagrangian is computed (has size *n* ).

new_x
*****
if true, no Ipopt evaluation method was previous called with the same
value for *x* .

obj_factor
**********
is the factor :math:`\alpha` that multiplies the objective f(x)
in the definition of the Lagrangian.

m
*
is the number of constraints in the problem (dimension of g(x)).

lambda
******
is the value of the constraint multipliers :math:`\lambda`
at which the Hessian is to be evaluated (has size *m* ).

new_lambda
**********
if true, no Ipopt evaluation method was previous called with the same
value for *lambda* .

nele_hess
*********
is the number of non-zero elements in the Hessian :math:`L_{x,x} (x)`; i.e.,
the same as
:ref:`ipopt_random_get_nlp_info@nnz_h_lag` .

iRow
****
If *values* is ``NULL`` ,
*iRow* has size *nele_hess* and is set to the
row indices for the non-zero entries in the
lower triangle of the Hessian :math:`L_{x,x} (x)`.

jCol
****
If *values* is ``NULL`` ,
*jCol* has size *nele_hess* and is set to the
column indices for the non-zero entries in the
lower triangle of the Hessian :math:`L_{x,x} (x)`.

values
******
If *values* is not ``NULL`` ,
it has size *nele_hess* and *values* [ *k* ]
is set to the value of element of the Hessian :math:`L_{x,x} (x)`
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

mixed_object\_.ran_hes_fun\_
****************************
After this call, the zero order Taylor coefficients in this function
will corresponding to the value of *fixed_vec_* and
the random effects in *x* .

Prototype
*********
{xrst_spell_off}
{xrst_code cpp} */
bool ipopt_random::eval_h(
   Index         n              ,  // in
   const Number* x              ,  // in
   bool          new_x          ,  // in
   Number        obj_factor     ,  // in
   Index         m              ,  // in
   const Number* lambda         ,  // in
   bool          new_lambda     ,  // in
   Index         nele_hess      ,  // in
   Index*        iRow           ,  // out
   Index*        jCol           ,  // out
   Number*       values         )  // out
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_random_eval_h}
*/
{  try
   {  try_eval_h(
         n,
         x,
         new_x,
         obj_factor,
         m,
         lambda,
         new_lambda,
         nele_hess,
         iRow,
         jCol,
         values
      );
   }
   catch(const std::exception& e)
   {  error_message_ = "ipopt_random::eval_h: std::exception: ";
      for(size_t j = 0; j < n_random_; j++)
         error_random_[j] = x[j];
      return false;
   }
   catch(const CppAD::mixed::exception& e)
   {  error_message_ = e.message("ipopt_random::eval_h");
      for(size_t j = 0; j < n_random_; j++)
         error_random_[j] = x[j];
      return false;
   }
   return true;
}
void ipopt_random::try_eval_h(
   Index         n              ,  // in
   const Number* x              ,  // in
   bool          new_x          ,  // in
   Number        obj_factor     ,  // in
   Index         m              ,  // in
   const Number* lambda         ,  // in
   bool          new_lambda     ,  // in
   Index         nele_hess      ,  // in
   Index*        iRow           ,  // out
   Index*        jCol           ,  // out
   Number*       values         )  // out
{  assert( size_t(n) == n_random_ );
   assert( m == 0 );
   assert( size_t(nele_hess) == nnz_h_lag_ );
   //
   if( new_x )
   {  // set the zero order Taylor coefficients in
      // mixed_object_.ran_like_fun_
      Number obj_value;
      eval_f(n, x, new_x, obj_value);
   }
   //
   const s_vector& row( mixed_object_.ran_hes_uu_rcv_.row() );
   const s_vector& col( mixed_object_.ran_hes_uu_rcv_.col() );
   assert( row.size() == nnz_h_lag_ );
   assert( col.size() == nnz_h_lag_ );
   //
   if( values == NULL )
   {  for(size_t k = 0; k < nnz_h_lag_; k++)
      {  // only returning lower triagle of Hessian of objective
         assert( col[k] <= row[k] );
         assert( row[k] < n_random_ );
         iRow[k] = static_cast<Index>( row[k] );
         jCol[k] = static_cast<Index>( col[k] );
      }
      assert( ! new_x );
      return;
   }
   //
   // random effects as a vector
   d_vector random_vec(n_random_);
   for(size_t j = 0; j < n_random_; j++)
      random_vec[j] = x[j];
   //
   // pack both the fixed and random effects into one vector
   d_vector both_vec(n_fixed_ + n_random_);
   mixed_object_.pack(fixed_vec_, random_vec, both_vec);
   //
   if( new_x )
   {  // set zero order Taylor coefficient in ran_like_fun_.
      d_vector vec = mixed_object_.ran_like_fun_.Forward(0, both_vec);
      if( CppAD::hasnan( vec ) ) throw CppAD::mixed::exception(
         "", "Hessian has a nan"
      );
   }
   //
   // computes the Hessian of objecive w.r.t random effects  f_uu (theta, u)
   d_vector val = mixed_object_.ran_hes_fun_.Forward(0, both_vec);
   if( CppAD::hasnan( val ) ) throw CppAD::mixed::exception(
      "", "Hessian has a nan"
   );
   assert( val.size() == nnz_h_lag_ );
   //
   // return the values
   for(size_t k = 0; k < nnz_h_lag_; k++)
      values[k] = obj_factor * static_cast<Number>( val[k] );
   //
   return;
}
} } // END_CPPAD_MIXED_NAMESPACE
