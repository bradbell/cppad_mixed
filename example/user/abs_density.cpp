// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin abs_density.cpp}

Absolute Value In Log-Density: Example and Test
###############################################

Model
*****

.. math::

   \B{p}( z_i | \theta ) \sim \B{L} ( \theta_i , \sigma )

where :math:`\B{L} ( \mu , \sigma )` is the Laplace distribution
with mean :math:`\mu` and standard deviation :math:`\sigma`.
The corresponding fixed likelihood
:ref:`g(theta)<theory@Fixed Likelihood, g(theta)>`
is

.. math::

   g( \theta ) = \sum_{i} \left[
      \log ( \sigma \sqrt{2} )
      +
      \sqrt{2} \; \left| \frac{ z_i - \exp( \theta_i )}{\sigma} \right|
   \right]

The optimal solution, with no constraints and no prior on :math:`\theta` is

.. math::

   \hat{\theta}_i = \log( z_i )

{xrst_literal
   // BEGIN C++
   // END C++
}

{xrst_end abs_density.cpp}
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
      double                sigma_;
      const d_vector&       z_;
   public:
      // constructor
      mixed_derived(
         size_t                 n_fixed       ,
         size_t                 n_random      ,
         bool                   quasi_fixed   ,
         bool                   bool_sparsity ,
         double                 sigma         ,
         const d_vector&        z             ) :
         cppad_mixed(n_fixed, n_random, quasi_fixed, bool_sparsity) ,
         n_fixed_(n_fixed)                                          ,
         sigma_(sigma)                                              ,
         z_(z)
      {  assert(z.size() == n_fixed); }

      // template version of fix_likelihood; i.e., p(z|theta) * p(theta)
      a1_vector fix_likelihood(
         const a1_vector&         fixed_vec  ) override
      {
         // initialize log-density
         a1_vector vec(1 + n_fixed_);
         vec[0] = 0.0;

         // compute this factors once
         a1_double sqrt_2 =  CppAD::sqrt( 2.0  );

         for(size_t j = 0; j < n_fixed_; j++)
         {  // Data term
            a1_double res   = z_[j] - CppAD::exp( fixed_vec[j] );
            res        /=  sigma_ ;
            // the following term does not depend on fixed effects
            // vec[0]     += log(sigma_ * sqrt_2);
            vec[1 + j] += sqrt_2 * res;
         }
         return vec;
      }
   };
}

bool abs_density_xam(void)
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
   d_vector random_lower(n_random), random_upper(n_random);
   std::string random_ipopt_options = "";
   //
   // no constriants
   d_vector fix_constraint_lower(0), fix_constraint_upper(0);
   //
   d_vector z(n_fixed);
   for(size_t i = 0; i < n_fixed; i++)
      z[i] = double(i+3);

   // object that is derived from cppad_mixed
   bool quasi_fixed   = false;
   bool bool_sparsity = true;
   double sigma       = 1.0;
   mixed_derived mixed_object(
      n_fixed, n_random, quasi_fixed, bool_sparsity, sigma, z
   );
   mixed_object.initialize(fixed_in, random_in);

   // optimize the fixed effects using quasi-Newton method
   std::string fixed_ipopt_options =
      "Integer print_level               0\n"
      "String  sb                        yes\n"
      "String  derivative_test           adaptive\n"
      "String  derivative_test_print_all yes\n"
      "Numeric tol                       1e-8\n"
      "Integer max_iter                  15\n"
   ;
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

   for(size_t j = 0; j < n_fixed; j++)
      ok &= fabs( fixed_out[j] - CppAD::log( z[j] ) ) <= tol;

   return ok;
}
// END C++
