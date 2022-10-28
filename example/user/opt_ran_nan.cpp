// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
$begin opt_ran_nan.cpp$$
$spell
   CppAD
   cppad
$$

$section Nan's During Optimization of Random Effects: Example and Test$$

$code
$srcthisfile%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>

namespace {
   using CppAD::vector;
   using CppAD::log;
   using CppAD::AD;
   //
   using CppAD::mixed::d_sparse_rcv;
   using CppAD::mixed::d_vector;
   //
   class mixed_derived : public cppad_mixed {
   private:
      const d_vector&       y_;
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
         ),
         y_(y)
      { }
      // implementation of ran_likelihood
      a1_vector ran_likelihood(
         const a1_vector&         theta  ,
         const a1_vector&         u      ) override
      {
         a1_vector vec(1);

         // initialize part of log-density that is always smooth
         vec[0] = 0.0;

         // sqrt_2pi = CppAD::sqrt( 8.0 * CppAD::atan(1.0) );

         // sum of residual squared
         a1_double sum_sq  = 0.0;
         for(size_t i = 0; i < y_.size(); i++)
         {  a1_double mu     = u[i];
            a1_double sigma  = theta[i];
            a1_double res    = (y_[i] - mu) / sigma;

            // Gaussian likelihood
            vec[0]  += log(sigma) + res * res / 2.0;
            // following term does not depend on fixed or random effects
            // vec[0]  += log(sqrt_2pi);

            // add to sum
            sum_sq += res * res;
         }

         // return nan when sum of squares is less than 1e-4
         a1_double vec_0 = CppAD::numeric_limits<a1_double>::quiet_NaN();
         a1_double small =  1e-4 ;
         vec_0  = CppAD::CondExpGt(sum_sq, + small, vec[0], vec_0);
         vec[0] = vec_0;
         //
         return vec;
      }
      // we expect to get a warnings
      virtual void warning(const std::string& warning_message)
      { }
   };
}

bool opt_ran_nan_xam(void)
{
   bool   ok = true;

   size_t n_data = 10;
   d_vector data(n_data), fixed_vec(n_data), random_in(n_data);

   for(size_t i = 0; i < n_data; i++)
   {  data[i]      = double(i + 1);
      fixed_vec[i] = 1.0;
      random_in[i] = 0.0;
   }

   // object that is derived from cppad_mixed
   bool quasi_fixed   = true;
   bool bool_sparsity = true;
   mixed_derived mixed_object(
      n_data, n_data, quasi_fixed, bool_sparsity, data
   );
   mixed_object.initialize(fixed_vec, random_in);

   // lower and upper limits for random effects
   double inf = std::numeric_limits<double>::infinity();
   d_vector random_lower(n_data), random_upper(n_data);
   for(size_t i = 0; i < n_data; i++)
   {  random_lower[i] = -inf;
      random_upper[i] = +inf;
   }

   // -----------------------------------------------------------------------
   // attempt to use ipopt to determine the optimal random effects
   std::string ipopt_options;
   ipopt_options += "Integer print_level      0\n";
   ipopt_options += "Integer max_iter         10\n";
   ipopt_options += "String  sb               yes\n";
   ipopt_options += "String  derivative_test  second-order\n";
   d_vector random_out = mixed_object.optimize_random(
      ipopt_options, fixed_vec, random_lower, random_upper, random_in
   );

   // check that the optimize had backed up and solve problem
   // to be near forbidden region (where nans occur)
   double sum_sq = 0.0;
   for(size_t i = 0; i < n_data; i++)
   {  double mu    = random_out[i];
      double sigma = fixed_vec[i];
      double res   = (data[i] - mu) / sigma;
      sum_sq      += res * res;
   }
   ok &= sum_sq >= 1e-4;
   ok &= sum_sq <= 1e-3;
   //
   return ok;
}
// END C++
