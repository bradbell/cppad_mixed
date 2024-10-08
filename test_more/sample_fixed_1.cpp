// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-24 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
$section Sample From Fixed Effects Posterior: Example and Test$$
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/manage_gsl_rng.hpp>
# include <Eigen/Dense>

namespace {
   using CppAD::vector;
   using CppAD::log;
   using CppAD::AD;
   //
   using CppAD::mixed::a1_double;
   using CppAD::mixed::d_sparse_rcv;
   //
   typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > double_mat;
   //
   // Used for debugging
   // void print(const char* name , const double_mat& mat)
   // {  std::cout << "\n" << name << " =\n" << mat << "\n"; }
   //
   class mixed_derived : public cppad_mixed {
   private:
      const vector<double>& y_;
   public:
      // constructor
      mixed_derived(
         size_t n_fixed                    ,
         size_t n_random                   ,
         bool   quasi_fixed                ,
         bool   bool_sparsity              ,
         const CppAD::mixed::d_sparse_rcv&    A_rcv  ,
         const vector<double>& y           ) :
         cppad_mixed(n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv),
         y_(y)
      {  assert( n_fixed      == y_.size() );
         assert( n_random % 2 == 0 );
         assert( n_random / 2 == y_.size() );
      }
      // implementation of ran_likelihood
      template <typename Vector>
      Vector template_ran_likelihood(
         const Vector& theta  ,
         const Vector& u      )
      {  typedef typename Vector::value_type scalar;

         Vector vec(1);

         // initialize part of log-density that is always smooth
         vec[0] = scalar(0.0);

         // pi
         scalar sqrt_2pi = scalar(
             CppAD::sqrt(8.0 * CppAD::atan(1.0)
         ));

         size_t n_data = y_.size();
         for(size_t i = 0; i < n_data; i++)
         {  scalar mu     = u[2*i] + u[2*i+1] + theta[i] + theta[0];
            scalar sigma  = scalar(1.0);
            scalar res    = (y_[i] - mu) / sigma;

            // p(y_i | u, theta)
            vec[0] += log(sqrt_2pi * sigma) + res*res / scalar(2.0);

            // p(u_i | theta)
            vec[0] += log(sqrt_2pi) + u[2*i] * u[2*i] / scalar(2.0);
            vec[0] += log(sqrt_2pi) + u[2*i+1] * u[2*i+1] / scalar(2.0);
         }
         return vec;
      }
      // a1_vector version of ran_likelihood
      virtual a1_vector ran_likelihood(
         const a1_vector& fixed_vec, const a1_vector& random_vec
      )
      {  return template_ran_likelihood( fixed_vec, random_vec ); }
      // ------------------------------------------------------------------
      // fix_constraint
      template <class Float>
      vector<Float> template_fix_constraint(
         const vector<Float>& fixed_vec  )
      {  vector<Float> ret_val(2);
         //
         for(size_t i = 0; i < 2; i++)
         {  ret_val[i] = 0.0;
            for(size_t j = 0; j < fixed_vec.size(); j++)
            {  if( i == 0 )
                  ret_val[i] += fixed_vec[j] * fixed_vec[j];
               else
                  ret_val[i] += CppAD::cos( fixed_vec[j] );
            }
         }
         return ret_val;
      }
   public:
      // ------------------------------------------------------------------
      // ran_likelihood
      // ------------------------------------------------------------------
      // fix_constraint
      a1_vector fix_constraint(
         const a1_vector& fixed_vec  )
      {  return template_fix_constraint(fixed_vec); }
      vector<double> fix_constraint(
         const vector<double>& fixed_vec  )
      {  return template_fix_constraint(fixed_vec); }
      // ------------------------------------------------------------------
   };
}

bool sample_fixed_1(void)
{
   bool   ok = true;
   double inf = std::numeric_limits<double>::infinity();
   //
   // initialize gsl random number generator
   size_t random_seed = CppAD::mixed::new_gsl_rng(0);
   //
   size_t n_data   = 5;
   size_t n_fixed  = n_data;
   size_t n_random = 2 * n_data;
   vector<double>
      fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
   for(size_t i = 0; i < n_fixed; i++)
   {  fixed_lower[i] = - inf;
      fixed_in[i]    = 0.1;
      fixed_upper[i] = inf;
   }
   //
   // explicit constriants (in addition to l1 terms)
   vector<double> fix_constraint_lower(2), fix_constraint_upper(2);
   for(size_t i = 0; i < 2; i++)
   {  fix_constraint_lower[i] = -inf;
      fix_constraint_upper[i] = inf;
   }
   //
   vector<double> data(n_data), random_in(n_random);
   for(size_t i = 0; i < n_data; i++)
   {  data[i]       = double(i + 1);
      random_in[i] = 0.0;
   }
   // constraint the sum of the first two random effect estimates to zero
   // (this results in a constraint on the first fixed effect)
   CppAD::mixed::sparse_rc A_pattern(1, n_random, 2);
   for(size_t k = 0; k < A_pattern.nnz(); k++)
      A_pattern.set(k, 0, k);
   CppAD::mixed::d_sparse_rcv A_rcv(A_pattern);
   for(size_t k = 0; k < A_rcv.nnz(); k++)
      A_rcv.set(k, 1.0);

   // object that is derived from cppad_mixed
   bool quasi_fixed = true;
   bool bool_sparsity = false;
   mixed_derived mixed_object(
      n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, data
   );
   mixed_object.initialize(fixed_in, random_in);

   // options for optimization
   std::string fixed_ipopt_options =
      "Integer print_level               0\n"
      "String  sb                        yes\n"
      "String  derivative_test           first-order\n"
      "String  derivative_test_print_all yes\n"
      "Numeric tol                       1e-8\n"
      "Integer max_iter                  15\n"
   ;
   std::string random_ipopt_options =
      "Integer print_level     0\n"
      "String  sb              yes\n"
      "String  derivative_test second-order\n"
      "Numeric tol             1e-8\n"
   ;
   // lower and upper limits for random effects
   vector<double> random_lower(n_random), random_upper(n_random);
   for(size_t i = 0; i < n_random; i++)
   {  random_lower[i] = -inf;
      random_upper[i] = +inf;
   }
   // optimize fixed effects
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
   //
   // now change the limit on the second fixed effect so it is active
   size_t bnd_j   = 1;
   fixed_lower[bnd_j] = 1.5 * solution.fixed_opt[bnd_j];
   fixed_in[bnd_j]    = 2.0 * solution.fixed_opt[bnd_j];
   //
   // also change the limit on the first fix_constraint so it is active
   vector<double> fix_con = mixed_object.fix_constraint(solution.fixed_opt);
   size_t con_i = 0;
   fix_constraint_lower[con_i] = 1.5 * fix_con[con_i];
   //
   // optimize fixed effects
   solution = mixed_object.optimize_fixed(
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
   //
   // check that one of the constraints is active
   ok &= solution.fixed_lag.size() == n_fixed;
   for(size_t j = 0; j < n_fixed; j++)
   {  if( j == bnd_j )
         ok &= solution.fixed_lag[j] != 0.0;
      else
         ok &= solution.fixed_lag[j] == 0.0;
   }
   ok &= solution.fix_con_lag.size() == 2;
   for(size_t i = 0; i < 2; i++)
   {  if( i == con_i )
         ok &= solution.fix_con_lag[i] != 0.0;
      else
         ok &= solution.fix_con_lag[i] == 0.0;
   }
   ok &= solution.ran_con_lag.size() == 1;
   ok &= solution.ran_con_lag[0] != 0.0;
   //
   // corresponding optimal random effects
   vector<double> random_opt = mixed_object.optimize_random(
      random_ipopt_options,
      solution.fixed_opt,
      random_lower,
      random_upper,
      random_in
   );
   //
   // compute corresponding information matrix
   d_sparse_rcv
   information_rcv = mixed_object.information_mat(solution, random_opt);
   //
   // sample from the posterior for fixed effects
   size_t n_sample = 10000;
   CppAD::vector<double> sample( n_sample * n_fixed );
   mixed_object.sample_fixed(
      sample,
      information_rcv,
      solution,
      fixed_lower,
      fixed_upper
   );
   //
   // compute sample covariance matrix
   double_mat sample_cov = double_mat::Zero(n_fixed, n_fixed);
   for(size_t i = 0; i < n_sample; i++)
   {  double_mat diff(n_fixed, 1);
      for(size_t j = 0; j < n_fixed; j++)
         diff(j, 0) = sample[ i * n_fixed + j] - solution.fixed_opt[j];
      sample_cov += diff * diff.transpose();
   }
   sample_cov *= 1.0 / double(n_sample);
   //
   // unconstrained information matrix
   double_mat info_mat = double_mat::Zero(n_fixed, n_fixed);
   size_t K = information_rcv.nnz();
   for(size_t k = 0; k < K; k++)
   {  size_t i = information_rcv.row()[k];
      size_t j = information_rcv.col()[k];
      info_mat(i, j) = information_rcv.val()[k];
      info_mat(j, i) = information_rcv.val()[k];
   }
   // unconstrained covraince
   double_mat C = info_mat.inverse();
   //
   double max_cov = 0.0;
   for(size_t i = 0; i < n_fixed; i++)
      max_cov = std::max( max_cov, C(i, i) );
   //
   double max_err = 0.0;
   for(size_t i = 0; i < n_fixed; i++)
   {  for(size_t j = 0; j < n_fixed; j++)
      {  double value = sample_cov(i, j);
         double check = C(i, j);
         max_err = std::max(max_err, std::fabs(value - check) / max_cov );
      }
   }
   ok &= max_err < .05;
   //
   if( ! ok )
   {  std::cout << "\nmax_err = " << max_err << "\n";
      std::cout << "\nrandom_seed = " << random_seed << "\n";
   }
   //
   CppAD::mixed::free_gsl_rng();
   return ok;
}
// END C++
