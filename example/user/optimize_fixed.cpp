// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin optimize_fixed.cpp}

Optimize Fixed Effects: Example and Test
########################################

Model
*****

.. math::

   \B{p}( y_i | \theta , u ) \sim \B{N} ( u_i + \theta_0 , \theta_1^2 )

.. math::

   \B{p}( u_i | \theta ) \sim \B{N} ( 0 , 1 )

.. math::

   \B{p}( \theta ) \sim \B{N} ( 4 , 1 )

It follows that the Laplace approximation is exact and

.. math::

   \B{p}( y_i | \theta ) \sim \B{N} \left( \theta_0 , 1 + \theta_1^2 \right)

The constraints on the fixed effect are

.. math::

   - \infty \leq \theta_0 \leq + \infty
   \R{\; and \;}
   0.1 \leq \theta_1 \leq 100

Objective
*********
The corresponding objective for the fixed effects is equivalent to:

.. math::

   F( \theta ) = \frac{1}{2} \left[
      ( \theta_0 - 4 )^2 + ( \theta_1 - 4 )^2 +
         N \log \left( 1 + \theta_1^2 \right) +
         ( 1 + \theta_1^2)^{-1} \sum_{i=0}^{N-1} ( y_i - \theta_0 )^2
   \right]

First Order Partials
********************
The first order partial derivatives of the objective are:

.. math::

   F_0 ( \theta )
   =
   ( \theta_0 - 4 )
   -
   ( 1 + \theta_1^2)^{-1} \sum_{i=0}^{N-1} ( y_i - \theta_0 )

.. math::

   F_1 ( \theta )
   =
   ( \theta_1 - 4 )
   +
   N \left( 1 + \theta_1^2 \right)^{-1} \theta_1
   -
   ( 1 + \theta_1^2)^{-2} \theta_1 \sum_{i=0}^{N-1} ( y_i - \theta_0 )^2

Optimizer Trace
***************
This example uses the optimizer trace information; see
:ref:`fixed_solution@trace_vec` .

Optimizer Warm Start
********************
This example uses the optimizer warm start information; see
:ref:`fixed_solution@warm_start` .

Source Code
***********
{xrst_literal
   // BEGIN C++
   // END C++
}

{xrst_end optimize_fixed.cpp}
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>

namespace {
   using CppAD::log;
   using CppAD::AD;
   //
   using CppAD::mixed::d_sparse_rcv;
   using CppAD::mixed::a1_double;
   using CppAD::mixed::d_vector;
   using CppAD::mixed::a1_vector;
   //
   class mixed_derived : public cppad_mixed {
   private:
      const size_t          n_fixed_;
      const size_t          n_random_;
      const d_vector&       y_;
   // ----------------------------------------------------------------------
   public:
      // constructor
      mixed_derived(
         size_t                 n_fixed       ,
         size_t                 n_random      ,
         bool                   quasi_fixed   ,
         bool                   bool_sparsity ,
         const d_vector&       y              ) :
         cppad_mixed(
            n_fixed, n_random, quasi_fixed, bool_sparsity
         )                     ,
         n_fixed_(n_fixed)     ,
         n_random_(n_random)   ,
         y_(y)
      {  assert( n_fixed == 2);
         assert( y_.size() == n_random_ );
      }
   // ----------------------------------------------------------------------
      // implementation of ran_likelihood
      a1_vector ran_likelihood(
         const a1_vector&         theta  ,
         const a1_vector&         u      ) override
      {
         assert( theta.size() == n_fixed_ );
         assert( u.size() == y_.size() );
         a1_vector vec(1);

         // initialize part of log-density that is always smooth
         vec[0] = 0.0;

         // sqrt_2pi = CppAD::sqrt(8.0 * CppAD::atan(1.0) );

         for(size_t i = 0; i < n_random_; i++)
         {  a1_double mu     = u[i] + theta[0];
            a1_double sigma  = theta[1];
            a1_double res    = (y_[i] - mu) / sigma;

            // p(y_i | u, theta)
            vec[0] += log(sigma) + res * res / 2.0;
            // following term does not depend on fixed or random effects
            // vec[0] += log(sqrt_2pi);

            // p(u_i | theta)
            vec[0] += u[i] * u[i] / 2.0;
            // following term does not depend on fixed or random effects
            // vec[0] += log(sqrt_2pi);
         }
         return vec;
      }
      // implementation of fix_likelihood
      a1_vector fix_likelihood(
         const a1_vector&         fixed_vec  ) override
      {
         assert( fixed_vec.size() == n_fixed_ );
         a1_vector vec(1);

         // initialize part of log-density that is smooth
         vec[0] = 0.0;

         // sqrt_2pi = CppAD::sqrt( 8.0 * CppAD::atan(1.0) );

         for(size_t j = 0; j < n_fixed_; j++)
         {  a1_double mu     = 4.0;
            a1_double sigma  = 1.0;
            a1_double res    = (fixed_vec[j] - mu) / sigma;

            // This is a Gaussian term, so entire density is smooth
            vec[0]  += res * res / 2.0;
            // following term does not depend on fixed effects
            // vec[0]  += log(sqrt_2pi * sigma);
         }
         return vec;
      }
   };
   // derivative of objective
   d_vector objective_fixed(
      const d_vector&       data   ,
      const d_vector&       theta  )
   {  d_vector dF(2);
      //
      // compute partials of F
      double sum   = 0.0;
      double sumsq = 0.0;
      for(size_t i = 0; i < data.size(); i++)
      {  sum   += theta[0] - data[i];
         sumsq += (theta[0] - data[i]) * (theta[0] - data[i]);
      }
      double den = 1.0 + theta[1] * theta[1];
      dF[0]  = (theta[0] - 4.0) + sum / den;
      dF[1]  = theta[1] - 4.0;
      dF[1] += double(data.size()) * theta[1] / den;
      dF[1] -= sumsq * theta[1]  / (den * den);
      //
      return dF;
   }
}

bool optimize_fixed_xam(void)
{
   bool   ok = true;
   double inf = std::numeric_limits<double>::infinity();
   double tol = 1e-8;

   size_t n_data   = 10;
   size_t n_fixed  = 2;
   size_t n_random = n_data;
   d_vector
      fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
   fixed_lower[0] = - inf; fixed_in[0] = 2.0; fixed_upper[0] = inf;
   fixed_lower[1] = .01;   fixed_in[1] = 0.5; fixed_upper[1] = inf;
   //
   // explicit constriants (in addition to l1 terms)
   d_vector fix_constraint_lower(0), fix_constraint_upper(0);
   //
   d_vector data(n_data), random_in(n_random);
   for(size_t i = 0; i < n_data; i++)
   {  data[i]       = double(i + 1);
      random_in[i] = 0.0;
   }

   // object that is derived from cppad_mixed
   bool quasi_fixed   = true;
   bool bool_sparsity = true;
   mixed_derived mixed_object(
         n_fixed, n_random, quasi_fixed, bool_sparsity, data
   );
   mixed_object.initialize(fixed_in, random_in);

   // optimize the fixed effects using quasi-Newton method
   std::string fixed_ipopt_options =
      "Integer print_level               0\n"
      "String  sb                        yes\n"
      "String  derivative_test           adaptive\n"
      "String  derivative_test_print_all yes\n"
      "Integer max_iter                  15\n"
   ;
   std::string random_ipopt_options =
      "Integer print_level     0\n"
      "String  sb              yes\n"
      "String  derivative_test second-order\n"
      "Numeric tol             1e-8\n"
   ;
   //
   d_vector random_lower(n_random), random_upper(n_random);
   for(size_t i = 0; i < n_random; i++)
   {  random_lower[i] = -inf;
      random_upper[i] = +inf;
   }
   // ------------------------------------------------------------------
   // optimize with tolerance 1e-3
   std::string temp_string = fixed_ipopt_options + "Numeric tol 1e-3\n";
   d_vector fixed_scale = fixed_in;
   CppAD::mixed::fixed_solution solution = mixed_object.optimize_fixed(
      temp_string,
      random_ipopt_options,
      fixed_lower,
      fixed_upper,
      fix_constraint_lower,
      fix_constraint_upper,
      fixed_scale,
      fixed_in,
      random_lower,
      random_upper,
      random_in
   );
   d_vector fixed_out = solution.fixed_opt;
   // ------------------------------------------------------------------
   // continue optimization, from previous, with new tolerance of 1e-8
   temp_string = fixed_ipopt_options + "Numeric tol 1e-8\n";
   fixed_in    = fixed_out;
   solution    = mixed_object.optimize_fixed(
      temp_string,
      random_ipopt_options,
      fixed_lower,
      fixed_upper,
      fix_constraint_lower,
      fix_constraint_upper,
      fixed_scale,
      fixed_in,
      random_lower,
      random_upper,
      random_in,
      solution.warm_start
   );
   fixed_out = solution.fixed_opt;
   // ------------------------------------------------------------------

   // deriative of objective at fixed_in and fixed_out
   d_vector dF_scale = objective_fixed(data, fixed_scale);
   d_vector dF_out   = objective_fixed(data, fixed_out);

   // scaling for objective
   double scale = std::max(
      std::fabs( dF_scale[0] ), std::fabs( dF_scale[1] )
   );
   scale = 1.0 / scale;

   // Note that no constraints are active, (not even the l1 terms)
   // so the partials should be zero.
   ok &= fabs( scale * dF_out[0] ) <= 5. * tol;
   ok &= fabs( scale * dF_out[1] ) <= 5. * tol;

   return ok;
}
// END C++
