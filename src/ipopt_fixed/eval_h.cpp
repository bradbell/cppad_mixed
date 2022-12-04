// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/ipopt_fixed.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
{xrst_begin ipopt_fixed_eval_h}
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
:ref:`ipopt_fixed_get_nlp_info@nnz_h_lag` .

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
if set to false, the optimization will terminate with status set to
:ref:`ipopt_fixed_finalize_solution@status@USER_REQUESTED_STOP` .

Prototype
*********
{xrst_spell_off}
{xrst_code cpp} */
bool ipopt_fixed::eval_h(
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

{xrst_end ipopt_fixed_eval_h}
*/
{
   double obj_factor_scaled;
   d_vector lambda_scaled(m);
   if( values == nullptr )
   {  double nan = std::numeric_limits<double>::quiet_NaN();
      obj_factor_scaled = nan;
      for(Index i = 0; i < m; ++i)
         lambda_scaled[i] = nan;
   }
   else
   {  for(Index j = 0; j < n; ++j)
         x_tmp_[j] = scale_x_[j] * x[j];
      obj_factor_scaled = scale_f_ * obj_factor;
      for(Index i = 0; i < m; i++)
         lambda_scaled[i] = scale_g_[i] * lambda[i];
   }
   if( abort_on_eval_error_ )
   {  try_eval_h(
         n,
         x_tmp_,
         new_x,
         obj_factor_scaled,
         m,
         lambda_scaled.data(),
         new_lambda,
         nele_hess,
         iRow,
         jCol,
         values
      );
   }
   else
   {  try
      {  try_eval_h(
            n,
            x_tmp_,
            new_x,
            obj_factor_scaled,
            m,
            lambda_scaled.data(),
            new_lambda,
            nele_hess,
            iRow,
            jCol,
            values
         );
      }
      catch(const std::exception& e)
      {  error_message_ = "ipopt_fixed::eval_h: std::exception: ";
         for(size_t j = 0; j < n_fixed_; j++)
            error_fixed_[j] = x[j];
         return false;
      }
      catch(const CppAD::mixed::exception& e)
      {  error_message_ = e.message("ipopt_fixed::eval_h");
         for(size_t j = 0; j < n_fixed_; j++)
            error_fixed_[j] = x[j];
         return false;
      }
   }
   assert( size_t( nele_hess ) == lag_hes_row_.size() );
   assert( size_t( nele_hess ) == lag_hes_col_.size() );
   if( values != nullptr )
   {  for(Index k = 0; k < nele_hess; ++k)
      {  size_t i = lag_hes_row_[k];
         size_t j = lag_hes_col_[k];
         values[k] *= scale_x_[i] * scale_x_[j];
      }
   }
   return true;
}
void ipopt_fixed::try_eval_h(
   Index           n              ,  // in
   const d_vector& x            ,  // in
   bool            new_x          ,  // in
   Number          obj_factor     ,  // in
   Index           m              ,  // in
   const Number*   lambda         ,  // in
   bool            new_lambda     ,  // in
   Index           nele_hess      ,  // in
   Index*          iRow           ,  // out
   Index*          jCol           ,  // out
   Number*         values         )  // out
{
   assert( ! mixed_object_.quasi_fixed_ );
   assert( n > 0 );
   assert( size_t(n) == n_fixed_ + fix_likelihood_nabs_ );
   assert( m >= 0 );
   assert( size_t(m) == 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_ );
   assert( size_t(nele_hess) == nnz_h_lag_ );
   if( values == NULL )
   {  for(size_t k = 0; k < nnz_h_lag_; k++)
      {  iRow[k] = Index( lag_hes_row_[k] );
         jCol[k] = Index( lag_hes_col_[k] );
      }
      return;
   }
   //
   // fixed effects
   for(size_t j = 0; j < n_fixed_; j++)
      fixed_tmp_[j] = double( x[j] );
   //
   // initialize return value
   for(size_t k = 0; k < nnz_h_lag_; k++)
      values[k] = Number( 0.0 );
   //
   // random part of objective
   if( n_random_ > 0 )
   {
      // compute the optimal random effects corresponding to fixed effects
      if( new_x )
         new_random(fixed_tmp_);
      // compute Hessian of random part of objective w.r.t. fixed effects
      w_laplace_obj_tmp_[0] = obj_factor;
      //
      // include random constraints in this Hessian calculation
      size_t offset = 2 * fix_likelihood_nabs_ + n_fix_con_;
      for(size_t i = 0; i < n_ran_con_; i++)
      {  w_laplace_obj_tmp_[i+1] = lambda[offset + i];
      }
      //
      mixed_object_.laplace_obj_hes(
         fixed_tmp_,
         random_cur_,
         w_laplace_obj_tmp_,
         laplace_obj_hes_info_.row,
         laplace_obj_hes_info_.col,
         laplace_obj_hes_info_.val
      );
      for(size_t k = 0; k < laplace_obj_hes_info_.row.size(); k++)
      {  size_t index = laplace_obj_hes_2_lag_[k];
         assert( index < nnz_h_lag_ );
         values[index] += Number( laplace_obj_hes_info_.val[k] );
      }
   }
   //
   // Hessian of Lagrangian of weighted fixed likelihood
   w_fix_likelihood_tmp_[0] = obj_factor;
   //
   for(size_t j = 0; j < fix_likelihood_nabs_; j++)
   {
      w_fix_likelihood_tmp_[1 + j] = lambda[2*j + 1] - lambda[2*j];
   }
   s_vector fix_like_hes_row = mixed_object_.fix_like_hes_.subset.row();
   s_vector fix_like_hes_col = mixed_object_.fix_like_hes_.subset.col();
   d_vector fix_like_hes_val = mixed_object_.fix_like_hes_.subset.val();
   mixed_object_.fix_like_hes(
      fixed_tmp_,
      w_fix_likelihood_tmp_,
      fix_like_hes_row,
      fix_like_hes_col,
      fix_like_hes_val
   );
   for(size_t k = 0; k < fix_like_hes_row.size(); k++)
   {  size_t index = fix_like_hes_2_lag_[k];
      assert( index < nnz_h_lag_ );
      values[index] += Number( fix_like_hes_val[k] );
   }
   //
   // Hessian of Lagrangian of fixed constraints
   for(size_t j = 0; j < n_fix_con_; j++)
   {  size_t ell        = 2 * fix_likelihood_nabs_ + j;
      w_fix_con_tmp_[j] = lambda[ell];
   }
   s_vector fix_con_hes_row = mixed_object_.fix_con_hes_.subset.row();
   s_vector fix_con_hes_col = mixed_object_.fix_con_hes_.subset.col();
   d_vector fix_con_hes_val = mixed_object_.fix_con_hes_.subset.val();
   mixed_object_.fix_con_hes(
      fixed_tmp_,
      w_fix_con_tmp_,
      fix_con_hes_row,
      fix_con_hes_col,
      fix_con_hes_val
   );
   for(size_t k = 0; k < fix_con_hes_row.size(); k++)
   {  size_t index = fix_con_hes_2_lag_[k];
      assert( index < nnz_h_lag_ );
      values[index] += Number( fix_con_hes_val[k] );
   }
   //
   for(size_t ell = 0; ell < nnz_h_lag_; ell++)
   {  if( CppAD::isnan( values[ell] ) ) throw CppAD::mixed::exception(
         "", "Hessian of Lagragian has a nan"
      );
   }
   assert( size_t(nele_hess) == nnz_h_lag_ );
   return;
}
} } // END_CPPAD_MIXED_NAMESPACE
