// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>

namespace {
   using CppAD::vector;
   using CppAD::log;
   using CppAD::AD;
   using CppAD::mixed::d_sparse_rcv;
   //
   using CppAD::mixed::a1_double;

   class mixed_derived : public cppad_mixed {
   private:
      size_t                n_fixed_;
      const vector<double>& z_;
   public:
      // constructor
      mixed_derived(
         size_t n_fixed                    ,
         size_t n_random                   ,
         bool   quasi_fixed                ,
         bool   bool_sparsity              ,
         const CppAD::mixed::d_sparse_rcv&    A_rcv,
         const vector<double>& z           ) :
         cppad_mixed(n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv),
         n_fixed_(n_fixed)                                      ,
         z_(z)
      {  assert(z.size() == n_fixed); }
      // implementation of fix_likelihood as p(z|theta) * p(theta)
      template <typename Vector>
      Vector template_fix_likelihood(
         const Vector& fixed_vec  )
      {  typedef typename Vector::value_type scalar;


         // initialize log-density
         Vector vec(1);
         vec[0] = scalar(0.0);

         for(size_t j = 0; j < n_fixed_ - 1; j++)
         {  // case where partial w.r.t. theta_j does not exist
            // when theta_j == z_j
            scalar sum_var;
            sum_var  = 0.0;
            //
            double    sum_data = 0.0;
            for(size_t k = 0; k <= j; k++)
            {  sum_var  += fixed_vec[k];
               sum_data += z_[k];
            }
            scalar res = sum_var - sum_data;
            vec[0] += 1e+12 * res * res / 2.0;
         }
         scalar res = z_[n_fixed_- 1] - fixed_vec[n_fixed_ - 1];
         vec[0]    += 1e+20 * res * res / 2.0;
         //
         return vec;
      }
      // a1_vector version of fix_likelihood
      virtual a1_vector fix_likelihood(const a1_vector& fixed_vec)
      {  return template_fix_likelihood( fixed_vec ); }
   };
}

bool scale_one(void)
{
   bool   ok = true;
   double inf = std::numeric_limits<double>::infinity();
   double tol = 1e-10;
   //
   // data
   size_t n_fixed  = 5;
   vector<double> z(n_fixed);
   double scale = 1.0;
   for(size_t i = 0; i < n_fixed; i++)
      z[i] = scale * double(i+1);

   // fixed effects
   vector<double>
      fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
   for(size_t j = 0; j < n_fixed; j++)
   {  fixed_lower[j] = - inf;
      fixed_in[j]    = 0.0;
      fixed_upper[j] = inf;
   }
   fixed_lower[n_fixed - 1] = 0.0;
   fixed_upper[n_fixed - 1] = 0.0;
   //
   // no random effects
   size_t n_random = 0;
   vector<double> random_in(0);
   //
   // no constriants
   vector<double> fix_constraint_lower(0), fix_constraint_upper(0);

   // object that is derived from cppad_mixed
   bool quasi_fixed = true;
   bool bool_sparsity = false;
   CppAD::mixed::d_sparse_rcv A_rcv; // empty matrix
   mixed_derived mixed_object(
      n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, z
   );
   mixed_object.initialize(fixed_in, random_in);

   // optimize the fixed effects using quasi-Newton method
   std::string fixed_ipopt_options =
      "Integer print_level                0\n"
      "String  sb                         yes\n"
      "String  derivative_test            first-order\n"
      "String  derivative_test_print_all yes\n"
      "Numeric tol                        1e-10\n"
      "Integer max_iter                   20\n"
      "Integer limited_memory_max_history 10\n"
   ;
   std::string random_ipopt_options =
      "Integer print_level 0\n"
      "String  sb          yes\n"
      "String  derivative_test second-order\n"
   ;
   vector<double> random_lower(n_random), random_upper(n_random);
   for(size_t i = 0; i < n_random; i++)
   {  random_lower[i] = -inf;
      random_upper[i] = +inf;
   }
   vector<double> fixed_scale = fixed_in;
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
   vector<double> fixed_out = solution.fixed_opt;
   //
   for(size_t j = 0; j < n_fixed - 1; j++)
   {  double err = fabs( fixed_out[j] / z[j] - 1.0 );
      ok &= err <= 5. * tol;
   }
   //
   ok &= fixed_out[n_fixed - 1] == 0.0;
   //
   return ok;
}
// END C++
