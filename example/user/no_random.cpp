// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin no_random.cpp}

No Random Effects: Example and Test
###################################

Model
*****

.. math::

   \B{p}( z_i | \theta ) \sim \B{N} ( \theta_i , 1 )

.. math::

   \B{p}( \theta_i ) \sim \B{N} ( 0 , 1 )

The corresponding fixed likelihood
:ref:`g(theta)<theory@Fixed Likelihood, g(theta)>`
is

.. math::

   g( \theta ) = \frac{1}{2} \sum_{i} \left[
      \log ( 2 \pi ) + \theta_i^2
      +
      \log ( 2 \pi ) + ( z_i - \theta_i )^2
   \right]

The optimal solution (with no constraints) is

.. math::

   \hat{\theta}_i = z_i / 2

{xrst_literal
   // BEGIN C++
   // END C++
}

{xrst_end no_random.cpp}
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
      size_t                n_fixed_;
      const d_vector&       z_;
   public:
      // constructor
      mixed_derived(
         size_t                 n_fixed        ,
         size_t                 n_random       ,
         const d_vector&        z              ) :
         cppad_mixed(n_fixed, n_random)  ,
         n_fixed_(n_fixed)               ,
         z_(z)
      {  assert(z.size() == n_fixed); }
      // implementation of fix_likelihood as p(z|theta) * p(theta)
      a1_vector fix_likelihood(
         const a1_vector&         fixed_vec  ) override
      {
         // initialize log-density
         a1_vector vec(1);
         vec[0] = 0.0;

         // compute this factors once
         // sqrt_2pi = CppAD::sqrt( 8.0 * CppAD::atan(1.0) );

         for(size_t j = 0; j < n_fixed_; j++)
         {
            // Data term p(z|theta)
            a1_double res  = (z_[j] - fixed_vec[j]);
            vec[0]    += res * res / 2.0;
            // following term does not depend on fixed effects
            // vec[0]    += log(sqrt_2pi );

            // prior term p(theta)
            res     = fixed_vec[j];
            vec[0] += res * res / 2.0;
            // following term does not depend on fixed effects
            // vec[0]    += log(sqrt_2pi );
         }
         return vec;
      }
   };
}

bool no_random_xam(void)
{
   bool   ok = true;
   double inf = std::numeric_limits<double>::infinity();
   double tol = 1e-8;

   // fixed effects
   size_t n_fixed  = 3;
   d_vector
      fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
   for(size_t j = 0; j < n_fixed; j++)
   {  fixed_lower[j] = - inf;
      fixed_in[j]    = 0.0;
      fixed_upper[j] = inf;
   }
   //
   // no random effects
   size_t n_random = 0;
   d_vector random_in(0);
   //
   // no constriants
   d_vector fix_constraint_lower(0), fix_constraint_upper(0);
   //
   d_vector z(n_fixed);
   for(size_t i = 0; i < n_fixed; i++)
      z[i] = double(i+1);

   // object that is derived from cppad_mixed
   // (test full netwon method to make sure it works with no random effects).
   mixed_derived mixed_object(n_fixed, n_random, z);
   mixed_object.initialize(fixed_in, random_in);

   // optimize the fixed effects using quasi-Newton method
   std::string fixed_ipopt_options =
      "Integer print_level               0\n"
      "String  sb                        yes\n"
      "String  derivative_test           first-order\n"
      "String  derivative_test_print_all yes\n"
      "Numeric tol                       1e-8\n"
      "Integer max_iter                  15\n"
   ;
   std::string random_ipopt_options =
      "Integer print_level 0\n"
      "String  sb          yes\n"
      "String  derivative_test second-order\n"
   ;
   d_vector random_lower(n_random), random_upper(n_random);
   for(size_t i = 0; i < n_random; i++)
   {  random_lower[i] = -inf;
      random_upper[i] = +inf;
   }
   d_vector fixed_scale = fixed_in;
   CppAD::mixed::fixed_solution solution = mixed_object.optimize_fixed(
      fixed_ipopt_options,
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
   //
   for(size_t j = 0; j < n_fixed; j++)
      ok &= fabs( fixed_out[j] - z[j] / 2.0 ) <= tol;

   return ok;
}
// END C++
