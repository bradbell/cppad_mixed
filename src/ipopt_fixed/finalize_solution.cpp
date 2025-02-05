// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-24 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/ipopt_fixed.hpp>
# include <IpIpoptData.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
// --------------------------------------------------------------------
// only used by finalize_solution
bool ipopt_fixed::check_in_limits(
   double lower, double x, double upper, double tol
)
{  // scale
   double scale = 0.0;
   if( lower != nlp_lower_bound_inf_ )
      scale = std::max(scale, std::fabs(lower) );
   if( upper != nlp_upper_bound_inf_ )
      scale = std::max(scale, std::fabs(upper) );
   if( scale == 0.0 )
      scale = 1.0;
   //
   // flag
   bool flag = true;
   flag     &= x <= upper + tol * scale;
   flag     &= lower - tol * scale <= x;
   //
   return flag;
}
/*
{xrst_begin ipopt_fixed_finalize_solution dev}
{xrst_spell
  cq
  curr
  doesn
  doxygen
  ip
  namespace
  naninf
  solution solution
}

Get Solution Results
####################

Syntax
******

| ``finalize_solution`` (
| |tab| *status* , *n* , *x* , *z_L* , *z_U* , *m* , *g* , *lambda* , *obj_value* , *ip_data* , *ip_cq*
| )

solution\_
**********
This routine checks the solution values and sets the member variable

   ``CppAD::mixed fixed_solution solution_``

n
*
is the number of variables in the problem (dimension of x).

x
*
is the final value (best value found) for the primal variables
(has size *n* ).

z_L
***
is the final value for the *x* lower bound constraint multipliers
(has size *n* ).

z_U
***
is the final value for the *x* upper bound constraint multipliers
(has size *n* ).

m
*
is the number of constraints in the problem (dimension of g(x)).

lambda
******
is the final value for the g(x) constraint multipliers :math:`\lambda`.

obj_value
*********
is the value of the objective f(x) at the final *x* value.

ip_data
*******
Expert users see the doxygen documentation for
``Ipopt::IpoptData`` class:

ip_data.curr_mu()
=================
This function call returns an ipopt ``Number`` equal to the
value of the penalty parameter *mu* .

ip_cq
*****
Unspecified; i.e., not part of the Ipopt user API.

status
******
These status values are in the ``Ipopt`` namespace; e.g.,
``SUCCESS`` is short for ``Ipopt::SUCCESS`` :

SUCCESS
=======
Algorithm terminated successfully at a locally optimal point,
satisfying the convergence tolerances (can be specified by options).

MAXITER_EXCEEDED
================
Maximum number of iterations exceeded (can be specified by an option).

CPUTIME_EXCEEDED
================
Maximum number of CPU seconds exceeded (can be specified by an option).

STOP_AT_TINY_STEP
=================
Algorithm proceeds with very little progress.

STOP_AT_ACCEPTABLE_POINT
========================
Algorithm stopped at a point that was converged, not to desired
tolerances, but to acceptable tolerances (see the acceptable-... options).

LOCAL_INFEASIBILITY
===================
Algorithm converged to a point of local infeasibility. Problem may be
infeasible.

USER_REQUESTED_STOP
===================
A user call-back function returned false, i.e.,
the user code requested a premature termination of the optimization.

DIVERGING_ITERATES
==================
It seems that the iterates diverge.

RESTORATION_FAILURE
===================
Restoration phase failed, algorithm doesn't know how to proceed.

ERROR_IN_STEP_COMPUTATION
=========================
An unrecoverable error occurred while Ipopt tried to compute
the search direction.

INVALID_NUMBER_DETECTED
=======================
Algorithm received an invalid number (such as NaN or Inf) from
the NLP; see also option check_derivatives_for_naninf.

Prototype
*********
{xrst_spell_off}
{xrst_code cpp} */
void ipopt_fixed::finalize_solution(
   Ipopt::SolverReturn               status    ,  // in
   Index                             n         ,  // in
   const Number*                     x         ,  // in
   const Number*                     z_L       ,  // in
   const Number*                     z_U       ,  // in
   Index                             m         ,  // in
   const Number*                     g         ,  // in
   const Number*                     lambda    ,  // in
   Number                            obj_value ,  // in
   const Ipopt::IpoptData*           ip_data   ,  // in
   Ipopt::IpoptCalculatedQuantities* ip_cq     )  // in
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_fixed_finalize_solution}
*/
{  bool ok = true;
   //
   assert( n > 0 );
   assert( size_t(n) == n_fixed_ + fix_likelihood_nabs_ );
   assert( m >= 0 );
   assert( size_t(m) == 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_ );
   //
   // get bounds information
   d_vector x_lower(n), x_upper(n), g_lower(m), g_upper(m);
   get_bounds_info(
      n, x_lower.data(), x_upper.data(), m, g_lower.data(), g_upper.data()
   );
   //
   // g_of_x
   d_vector g_of_x(m);
   bool new_x = true;
   eval_g(n, x, new_x, m, g_of_x.data());
   new_x = false;
   //
   // solution_.fixed_opt
   assert( solution_.fixed_opt.size() == 0 );
   solution_.fixed_opt.resize(n_fixed_);
   for(size_t j = 0; j < n_fixed_; j++)
      solution_.fixed_opt[j] = scale_x_[j] * x[j];
   //
   // short name for fixed effects tolerance
   double tol = fixed_tolerance_;
   //
   // check that x is within its limits
   for(size_t j = 0; j < size_t(n); j++)
   {  ok &= check_in_limits(
         x_lower[j], x[j], x_upper[j], 2.0 * tol
      );
   }
   //
   // check that g is within its limits
   for(size_t i = 0; i < size_t(m); i++)
   {  ok &= check_in_limits(
         g_lower[i], g[i], g_upper[i], 2.0 * tol
      );
   }
   //
   // check that the bound multipliers are feasible
   for(size_t j = 0; j < size_t(n); j++)
   {  ok &= 0.0 <= z_L[j];
      ok &= 0.0 <= z_U[j];
   }
   //
   // grad_f
   CppAD::vector<Number> grad_f(n);
   eval_grad_f(n, x, new_x, grad_f.data() );

   // jac_g
   CppAD::vector<Number> jac_g(nnz_jac_g_);
   eval_jac_g(
      n,
      x,
      new_x,
      m,
      Index(nnz_jac_g_),
      nullptr,
      nullptr,
      jac_g.data()
   );
   //
   // solution_.ran_con_lag
   assert( solution_.ran_con_lag.size() == 0 );
   solution_.ran_con_lag.resize(n_ran_con_);
   size_t offset = 2 * fix_likelihood_nabs_ + n_fix_con_;
   for(size_t j = 0; j < n_ran_con_; j++)
   {  assert( g_lower[offset + j] == g_upper[offset + j] );
      solution_.ran_con_lag[j] = lambda[ offset + j];
   }
   //
   // check complementarity conditions for gl <= g(x) <= gu
   // set solution_.fix_con_lag
   solution_.fix_con_lag.resize(n_fix_con_);
   for(size_t i = 0; i < size_t(m); ++i)
   {  double abs_lam_i = std::fabs( lambda[i] );
      //
      bool at_lower = false;
      if( g_lower[i] != nlp_lower_bound_inf_ && lambda[i] < 0.0 )
      {  ok &= (g[i] - g_lower[i]) * abs_lam_i < 10.0 * tol;
         at_lower = g[i] - g_lower[i] < abs_lam_i;
      }
      //
      bool at_upper = false;
      if( g_upper[i] != nlp_upper_bound_inf_ && lambda[i] > 0.0)
      {  ok &= (g_upper[i] - g[i]) * abs_lam_i < 10.0 * tol;
         at_upper = g_upper[i] - g[i] < abs_lam_i;
      }
      offset = 2 * fix_likelihood_nabs_;
      if( offset <= i && i < offset + n_fix_con_ )
      {  size_t j = i - offset;
         if( ! (at_lower || at_upper) )
            solution_.fix_con_lag[j] = 0.0;
         else
            solution_.fix_con_lag[j] = lambda[i];
      }
   }

   // check gradient of Lagragian L(x)
   // solution_.fixed_lag
   assert( solution_.fixed_lag.size() == 0 );
   solution_.fixed_lag.resize(n_fixed_);
   for(size_t j = 0; j < size_t(n); ++j)
   {  // sum will accumulate partial of L(x) w.r.t. x[j]
      double sum = grad_f[j];
      for(size_t k = 0; k < nnz_jac_g_; k++)
      {  if( jac_g_col_[k] == j )
         {  size_t i = jac_g_row_[k];
            sum     += lambda[i] * jac_g[k];
         }
      }
      //
      // at_lower, ok
      bool at_lower = false;
      if( x_lower[j] != nlp_lower_bound_inf_ )
      {  ok      &= (x[j] - x_lower[j]) * z_L[j] < 10.0 * tol;
         at_lower = x[j] - x_lower[j] <= z_L[j];
      }
      //
      // at_upper, ok
      bool at_upper = false;
      if( x_upper[j] != nlp_upper_bound_inf_ )
      {  ok &= (x_upper[j] - x[j]) * z_U[j] < 10.0 * tol;
         at_lower = x_upper[j] - x[j] <= z_U[j];
      }
      //
      // solution_.fixed_lag[j]
      if( j < n_fixed_ )
      {  if( at_lower || at_upper )
            solution_.fixed_lag[j] = - sum;
         else
            solution_.fixed_lag[j] = 0.0;
      }
      //
      if( x_lower[j] == x_upper[j] )
      {   // this constriant gets removed
         sum = 0.0;
      }
      else
      {  sum -= z_L[j];
         sum += z_U[j];
      }
      //
      ok &= std::fabs(sum) < 10.0 * tol;
   }

   // set warm_start information
   solution_.warm_start.mu                    = ip_data->curr_mu();
   solution_.warm_start.scale_f               = scale_f_;
   solution_.warm_start.x_info.resize(n);
   for(size_t j = 0; j < size_t(n); ++j)
   {  solution_.warm_start.x_info[j].x       = x[j];
      solution_.warm_start.x_info[j].z_L     = z_L[j];
      solution_.warm_start.x_info[j].z_U     = z_U[j];
      solution_.warm_start.x_info[j].scale_x = scale_x_[j];
   }
   solution_.warm_start.g_info.resize(m);
   for(size_t i = 0; i < size_t(m); ++i)
   {  solution_.warm_start.g_info[i].lambda   = lambda[i];
      solution_.warm_start.g_info[i].scale_g  = scale_g_[i];
   }

   // set member variable finalize_solution_ok_
   finalize_solution_ok_ = ok;
}
} } // END_CPPAD_MIXED_NAMESPACE
