// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/ipopt_fixed.hpp>
# include <cppad/mixed/one_dim_derivative_chk.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

/*
{xrst_begin ipopt_fixed_set_scaling dev}

Set Scaling Factors
###################

Syntax
******

| *ok* = ``set_scaling`` (
| |tab| *x_scale* , *x_lower* , *x_upper* , *grad_f* , *jac_g*
| )

x_scale
*******
is the argument value at which the scaling is computer.

x_lower
*******
is the optimization lower limit for x.

x_upper
*******
is the optimization upper limit for x.

grad_f
******
is the gradient of f(x) at *x_scale* .

jac_g
*****
is the Jacobian of g(x) at *x_scale* .
This is a sparse representation of the jacobian using the member variables
``jac_g_row_`` and ``jac_g_col_`` .

scale_f\_
*********
The input value of this member variable must be one.
Upon return it is set to a value that should be used to multiply
the values return by ``eval_f`` .

scale_g\_
*********
The input value of this member variable must be a vector with size
equal to the range space of g(x) and all its components must be one.
Upon return ``scale_g_`` [ *i* ] is set to a value that should be used to
multiply the *i*-th component returned by ``eval_g`` .

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE
   // END_PROTOTYPE
}

{xrst_end ipopt_fixed_set_scaling}
-------------------------------------------------------------------------------
*/
// BEGIN_PROTOTYPE
bool ipopt_fixed::set_scaling(
   const d_vector& x_scale ,
   const d_vector& x_lower ,
   const d_vector& x_upper ,
   const d_vector& grad_f  ,
   const d_vector& jac_g   )
// END_PROTOTYPE
{  //
   // m: number of components if g
   size_t m    = 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_;
   //
# ifndef NDEBUG
   // n: number of arguemnts to function being optimized
   size_t  n   = n_fixed_ + fix_likelihood_nabs_;
   //
   assert( scale_f_ == 1.0 );
   assert( scale_g_.size() == m );
   assert( scale_x_.size() == n );
   for(size_t i = 0; i < m; ++i)
      assert( scale_g_[i] == 1.0 );
   for(size_t j = 0; j < n; ++j)
      assert( scale_x_[j] == 1.0 );
# endif
   //
   // scale_max: maximum scaling factor
   const double scale_max = 1e+14;
   //
   // scale_min: minimum scaling factor
   double scale_min       = 1.0 / scale_max;
   // --------------------------------------------------------------------
   // scale_x_: Note yet working
   // --------------------------------------------------------------------
   /*
   d_vector max_partial(n);
   for(size_t j = 0; j < n; ++j)
      max_partial[j] = std::fabs( grad_f[j] );
   for(size_t k = 0; k < nnz_jac_g_; k++)
   {  size_t j = jac_g_col_[k];
         max_partial[j] = std::max(max_partial[j], std::fabs( jac_g[k] ) );
   }
   for(size_t j = 0; j < n; ++j)
   {  scale_x_[j] = 1.0 / max_partial[j];
      scale_x_[j] = std::max(scale_x_[j], scale_min);
      scale_x_[j] = std::min(scale_x_[j], scale_max);
   }
   */
   // --------------------------------------------------------------------
   // norm_jac_g[i] : maximum absolute partial of constraint component g[i],
   // including scale_x and not counting components with equal limits.
   // --------------------------------------------------------------------
   d_vector norm_jac_g(m);
   for(size_t i = 0; i < m; i++)
      norm_jac_g[i] = 0.0;
   for(size_t k = 0; k < nnz_jac_g_; k++)
   {  size_t i = jac_g_row_[k];
      size_t j = jac_g_col_[k];
      if( x_lower[j] < x_upper[j] )
      {  double scaled  = jac_g[k] * scale_x_[j];
         norm_jac_g[i] = std::max( norm_jac_g[i], std::fabs( scaled ) );
      }
   }
   // -----------------------------------------------------------------------
   // norm_grad_f : maximum absolute partial of objective f(x),
   // including scale_x and not counting components with equal limits.
   // Use corresponding components in g(x) for absolute value terms.
   // -----------------------------------------------------------------------
   double norm_grad_f = 0.0;
   for(size_t j = 0; j < n_fixed_; j++) if( x_lower[j] < x_upper[j] )
   {  double scaled = grad_f[j] * scale_x_[j];
      norm_grad_f = std::max( norm_grad_f, std::fabs(scaled) );
   }
# ifndef NDEBUG
   // skip axuillary variable terms (gradients are one)
   for(size_t j = 0; j < fix_likelihood_nabs_; j++)
      assert( grad_f[n_fixed_ + j] == 1.0 );
# endif
   // include absolute value terms in g(x) in norm_grad_f
   for(size_t i = 0; i < fix_likelihood_nabs_; i++)
   {  assert( norm_jac_g[2*i] == norm_jac_g[2*i+1] );
      norm_grad_f = std::max( norm_grad_f, std::fabs( norm_jac_g[2*i] ) );
   }
   // ----------------------------------------------------------------------
   // scale_f_
   // ----------------------------------------------------------------------
   scale_f_ = 1.0 / norm_grad_f;
   scale_f_ = std::max( scale_min, scale_f_);
   scale_f_ = std::min( scale_max, scale_f_);
   // -----------------------------------------------------------------------
   // set scale_g_
   // -----------------------------------------------------------------------
   assert( scale_g_.size() == m );
   for(size_t i = 0; i < m; i++)
   {  scale_g_[i] = 1.0 / norm_jac_g[i];
      scale_g_[i] = std::max(scale_min, scale_g_[i] );
      scale_g_[i] = std::min(scale_max, scale_g_[i] );
   }
   //
   return true;
}
} } // END_CPPAD_MIXED_NAMESPACE
