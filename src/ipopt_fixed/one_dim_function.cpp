// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/ipopt_fixed.hpp>
# include <cppad/mixed/one_dim_derivative_chk.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

/*
{xrst_begin ipopt_fixed_one_dim_function}
{xrst_spell
   enum
   lagrange
}
Compute Model Functions Along One Direction
###########################################

Syntax
******

   *ok* = ``one_dim_function`` ( *x_j* , *fun_out* )

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE
   // END_PROTOTYPE
}

Member Variables
****************

one_dim_function_j\_
====================
This determines which component of the optimization variables
*x* we are evaluating this function along.
This is referred to as the index *j* below.

one_dim_function_x\_
====================
This determines the value for *x* , except for the *j*-th component,
that is used to evaluate the function.
It's *j*-th component is changed to *x_j* and than
restored back to it's original value before this function returns.

one_dim_function_eval\_
=======================
This determines which function is being evaluated.
Its possible values are:

.. csv-table::
   :widths: auto

   ``eval_f_enum``,objective function f(x)
   ``eval_g_enum``,constraint function g(x)
   ``eval_grad_L_enum``,gradient of the Lagrangian L'(x)

one_dim_function_obj_factor\_
=============================
This determines the objective function and is referred to
as :math:`\alpha` in the definition of L(x) below.

one_dim_function_lambda\_
=========================
This determines the Lagrange multipliers and is referred to
as :math:`\lambda` in the definition of L(x) below.

x_j
***
This is the value used for the *j*-th component of the optimization variables
when evaluating the function; see ``one_dim_function_j_`` above.

fun_out
*******
The function value is returned in this vector.
It has the following size depending on the value of
*one_dim_function_eval_* :

.. csv-table::
   :widths: auto

   *one_dim_function_eval_*,*fun_out* . ``size`` ()
   ``eval_f_enum``,size of range space for f(x); i.e. 1
   ``eval_g_enum``,size of range space for g(x); i.e. *m*
   ``eval_grad_L_enum``,number of optimization variables; i.e. *n*

The input value of its elements does not matter.
Upon return it contains the specified function value.

Definition of L(x)
******************

.. math::

   L(x) = \alpha f(x) + \sum_{i=1}^m \lambda_i g_i (x)

{xrst_end ipopt_fixed_one_dim_function}
*/
// BEGIN_PROTOTYPE
bool ipopt_fixed::one_dim_function(double x_j, d_vector& fun_out)
// END_PROTOTYPE
{  // new_x
   bool    new_x          = true;
   //
   // m: number of components if g
   size_t m    = 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_;
   //
   // n: number of arguemnts to function being optimized
   size_t  n   = n_fixed_ + fix_likelihood_nabs_;
   //
   // j: component in argument space were one dimensional function defined
   size_t  j   = one_dim_function_j_;
   //
   // x_save: original value of j-th component of x (used to restore value)
   double  x_save = one_dim_function_x_[j];
   //
   // replace the j-th compponent with requested value
   one_dim_function_x_[j] = x_j;
   //
   // x: pointer to argument for n dimensional function
   Number* x  = one_dim_function_x_.data();
   //
   // ok: did the function evaluation complete
   bool ok = false;
   //
   switch( one_dim_function_eval_ )
   {  // evaluting f(x) along j-th component of x
      case eval_f_enum:
      assert( fun_out.size() == 1 );
      ok = eval_f(Index(n), x, new_x, fun_out[0]);
      break;
      // --------------------------------------------------------------------
      // evaluating g(x) along j-th component of x
      case eval_g_enum:
      assert( fun_out.size() == m );
      ok = eval_g(
         Index(n), x, new_x, Index(m), fun_out.data()
      );
      break;
      // --------------------------------------------------------------------
      // evaluating grad L(x) along j-th component of x
      case eval_grad_L_enum:
      assert( fun_out.size() == n );
      //
      // grad_f
      d_vector grad_f(n);
      ok  = eval_grad_f(
         Index(n), x, new_x, grad_f.data()
      );
      //
      // jac_g
      new_x       = false;
      Index* iRow = nullptr;
      Index* jCol = nullptr;
      d_vector jac_g(nnz_jac_g_);
      ok &= eval_jac_g(
         Index(n),
         x,
         new_x,
         Index(m),
         Index( nnz_jac_g_ ),
         iRow,
         jCol,
         jac_g.data()
      );
      //
      // obj_factor:
      double obj_factor = one_dim_function_obj_factor_;
      //
      // lambda: Lagrange multipliers,`
      d_vector& lambda  = one_dim_function_lambda_;
      //
      // L(x) = obj_factor * f(x) = sum_i lambda[i] * g_i(x)
      //
      // fun_out = gradient of L(x)
      for(size_t j1 = 0; j1 < n; ++j1)
         fun_out[j1] = obj_factor * grad_f[j1];
      for(size_t k = 0; k < nnz_jac_g_; ++k)
      {  size_t i     = jac_g_row_[k];
         size_t j1    = jac_g_col_[k];
         fun_out[j1] += lambda[i] * jac_g[k];
      }
      break;
      // --------------------------------------------------------------------
   }
   // restore one_dimensional_function_x_
   one_dim_function_x_[j] = x_save;
   //
   // return function evaluation status
   return ok;
}

} } // END_CPPAD_MIXED_NAMESPACE
