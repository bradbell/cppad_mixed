// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/ipopt_fixed.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
// --------------------------------------------------------------------
/*
{xrst_begin ipopt_fixed_fixed_eq_constrain dev}

Set Solution When All Fixed Effects Are Equality Constrained
############################################################

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE
   // END_PROTOTYPE
}

fixed_opt
*********
This is equal to the lower and upper limits for the fixed effects.

solution\_
**********
This member variable must be empty when ``fixed_eq_constrain`` is called
(see assertions above).

solution\_.fixed_opt
********************
This vector is set equal to *fixed_opt*.

solution\_.ran_con_lag
**********************
This vector is set to have size n_ran_con\_ and value zero.

solution\_.fix_con_lag
**********************
This vector is set to have size n_fix_con\_ and value zero.

solution\_.fixed_lag
********************
This vector is set to have size n_fixed\_ and value equal to the
negative of the gradient of the objective function f(x).
This gradient is in ipopt_fixed scaled coordinates
and should probably be in fixed effects coordinates.
(A similar problem exists in :ref:`ipopt_fixed_finalize_solution-name` .)


{xrst_end ipopt_fixed_fixed_eq_constrain}
*/
// BEGIN_PROTOTYPE
void ipopt_fixed::fixed_eq_constrain(const d_vector& fixed_opt )
{  assert( fixed_opt.size() == n_fixed_ );
   assert( solution_.fixed_opt.size() == 0 );
   assert( solution_.ran_con_lag.size() == 0 );
   assert( solution_.fixed_lag.size() == 0 );
   assert( solution_.warm_start.x_info.size() == 0 );
   assert( solution_.warm_start.g_info.size() == 0 );
// END_PROTOTYPE
   //
   // solution_.fixed_opt
   solution_.fixed_opt.resize(n_fixed_);
   for(size_t j = 0; j < n_fixed_; ++j)
      solution_.fixed_opt[j] = fixed_opt[j];
   //
   // solution_.ran_con_lag
   solution_.ran_con_lag.resize(n_ran_con_);
   for(size_t j = 0; j < n_ran_con_; j++)
      solution_.ran_con_lag[j] = 0.0;

   // solution_.fix_con_lag
   solution_.fix_con_lag.resize(n_fix_con_);
   for(size_t i = 0; i < size_t(n_fix_con_); ++i)
      solution_.fix_con_lag[i] = 0.0;

   // grad_f
   d_vector grad_f(n_fixed_), x_opt(n_fixed_);
   bool new_x = true;
   for(size_t j = 0; j < n_fixed_; ++j)
      x_opt[j] = fixed_opt[j] / scale_x_[j];
   int n = int( n_fixed_ );
   eval_grad_f(n, x_opt.data(), new_x, grad_f.data() );

   // solution_.fixed_lag
   solution_.fixed_lag.resize(n_fixed_);
   for(size_t j = 0; j < n_fixed_; ++j)
      solution_.fixed_lag[j] = - grad_f[j];
}

// --------------------------------------------------------------------
} } // END_CPPAD_MIXED_NAMESPACE
