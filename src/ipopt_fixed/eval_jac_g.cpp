// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/ipopt_fixed.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
{xrst_begin ipopt_fixed_eval_jac_g dev}
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
:ref:`ipopt_fixed_get_nlp_info@nnz_jac_g` .

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
if set to false, the optimization will terminate with status set to
:ref:`ipopt_fixed_finalize_solution@status@USER_REQUESTED_STOP` .

Prototype
*********
{xrst_spell_off}
{xrst_code cpp} */
bool ipopt_fixed::eval_jac_g(
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

{xrst_end ipopt_fixed_eval_jac_g}
*/
{  if( values != nullptr )
   {  for(Index j = 0; j < n; ++j)
         x_tmp_[j] = scale_x_[j] * x[j];
   }
   if( abort_on_eval_error_ )
   {  try_eval_jac_g(n, x_tmp_, new_x, m, nele_jac, iRow, jCol, values);
   }
   else
   {  try
      {  try_eval_jac_g(n, x_tmp_, new_x, m, nele_jac, iRow, jCol, values);
      }
      catch(const std::exception& e)
      {  error_message_ = "ipopt_fixed::eval_jac_g: std::exception: ";
         for(size_t j = 0; j < n_fixed_; j++)
            error_fixed_[j] = x[j];
         return false;
      }
      catch(const CppAD::mixed::exception& e)
      {  error_message_ = e.message("ipopt_fixed::eval_jac_g");
         for(size_t j = 0; j < n_fixed_; j++)
            error_fixed_[j] = x[j];
         return false;
      }
   }
   if( values != nullptr )
   {  for(size_t ell = 0; ell < nnz_jac_g_; ell++)
      {  size_t i     = jac_g_row_[ell];
         size_t j     = jac_g_col_[ell];
         values[ell] *= scale_g_[i] * scale_x_[j];
      }
   }
   return true;
}
void ipopt_fixed::try_eval_jac_g(
   Index           n        ,  // in
   const d_vector& x        ,  // in
   bool            new_x    ,  // in
   Index           m        ,  // in
   Index           nele_jac ,  // in
   Index*          iRow     ,  // out
   Index*          jCol     ,  // out
   Number*         values   )  // out
{
   assert( n > 0 );
   assert( size_t(n) == n_fixed_ + fix_likelihood_nabs_ );
   assert( m >= 0 );
   assert( size_t(m) == 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_ );
   assert( size_t(nele_jac) == nnz_jac_g_ );
   //
   s_vector fix_like_jac_row = mixed_object_.fix_like_jac_.subset.row();
   s_vector fix_like_jac_col = mixed_object_.fix_like_jac_.subset.col();
   d_vector fix_like_jac_val = mixed_object_.fix_like_jac_.subset.val();
   s_vector fix_con_jac_row = mixed_object_.fix_con_jac_.subset.row();
   s_vector fix_con_jac_col = mixed_object_.fix_con_jac_.subset.col();
   d_vector fix_con_jac_val = mixed_object_.fix_con_jac_.subset.val();
   if( values == NULL )
   {  // just return row and column indices for l1 constraints
      size_t ell = 0;
      for(size_t k = 0; k < fix_like_jac_row.size(); k++)
      {  if( fix_like_jac_row[k] != 0 )
         {  assert( ell + 1 < nnz_jac_g_ );
            iRow[ell] = Index( 2 * fix_like_jac_row[k] - 2 );
            jCol[ell] = Index( fix_like_jac_col[k] );
            ell++;
            iRow[ell] = Index( 2 * fix_like_jac_row[k] - 1 );
            jCol[ell] = Index( fix_like_jac_col[k] );
            ell++;
         }
      }
      // auxillary variables for l1 constraints
      for(size_t j = 0; j < fix_likelihood_nabs_; j++)
      {  assert( ell + 1 < nnz_jac_g_ );
         iRow[ell] = Index( 2 * j );
         jCol[ell] = Index( n_fixed_ + j);
         ell++;
         iRow[ell] = Index( 2 * j + 1);
         jCol[ell] = Index(n_fixed_ + j);
         ell++;
      }
      // fixed constraints
      size_t offset = 2 * fix_likelihood_nabs_;
      for(size_t k = 0; k < fix_con_jac_row.size(); k++)
      {  assert( ell < nnz_jac_g_ );
         iRow[ell] = Index( offset + fix_con_jac_row[k] );
         jCol[ell] = Index( fix_con_jac_col[k] );
         ell++;
      }
      // random constraints
      offset = 2 * fix_likelihood_nabs_ + n_fix_con_;
      for(size_t k = 0; k < ran_con_jac_rcv_.nnz(); k++)
      {  assert( ell < nnz_jac_g_ );
         iRow[ell] = Index( offset + ran_con_jac_rcv_.row()[k] );
         jCol[ell] = Index( ran_con_jac_rcv_.col()[k] );
         ell++;
      }
      assert( ell == nnz_jac_g_ );
      //
      return;
   }
   assert( jac_g_row_.size() == size_t(nnz_jac_g_) );
   //
   // fixed effects
   for(size_t j = 0; j < n_fixed_; j++)
      fixed_tmp_[j] = double( x[j] );
   //
   // check if this is a new x
   if( new_x && n_random_ > 0 )
      new_random(fixed_tmp_);
   //
   // Jacobian of fixed effects likelihood
   // (2DO: do not revaluate when eval_grad_f had same x)
   mixed_object_.fix_like_jac(
      fixed_tmp_,
      fix_like_jac_row,
      fix_like_jac_col,
      fix_like_jac_val
   );
   size_t ell = 0;
   for(size_t k = 0; k < fix_like_jac_row.size(); k++)
   {  if( fix_like_jac_row[k] != 0 )
      {  assert( ell + 1 < nnz_jac_g_ );
         values[ell] = Number( - fix_like_jac_val[k] );
         ell++;
         values[ell] = Number( + fix_like_jac_val[k] );
         ell++;
      }
   }
   for(size_t j = 0; j < fix_likelihood_nabs_; j++)
   {  assert( ell + 1 < nnz_jac_g_ );
      values[ell+1] = values[ell] = Number(1.0);
      ell += 2;
   }
   //
   // Jacobian of fixed constraints
   mixed_object_.fix_con_jac(
      fixed_tmp_,
      fix_con_jac_row,
      fix_con_jac_col,
      fix_con_jac_val
   );
   for(size_t k = 0; k < fix_con_jac_row.size(); k++)
   {  assert( ell < nnz_jac_g_ );
      values[ell++] = Number( fix_con_jac_val[k] );
   }
   //
   // Jacobian of random constraints
   if( n_ran_con_ > 0 )
   {  assert( n_random_ > 0 );
      mixed_object_.ran_con_jac(fixed_tmp_, random_cur_, ran_con_jac_rcv_);
      for(size_t k = 0; k < ran_con_jac_rcv_.nnz(); k++)
      {  assert( ell < nnz_jac_g_ );
         values[ell++] = Number( ran_con_jac_rcv_.val()[k] );
      }
   }
   assert( ell == nnz_jac_g_ );
   //
   for(ell = 0; ell < nnz_jac_g_; ell++)
   {  if( CppAD::isnan( values[ell] ) ) throw CppAD::mixed::exception(
         "", "constraint Jacobian has a nan"
      );
   }
   return;
}
} } // END_CPPAD_MIXED_NAMESPACE
