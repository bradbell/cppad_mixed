// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/sparse_mat_info.hpp>
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
   class mixed_derived : public cppad_mixed {
   private:
      const size_t          n_fixed_;
      const size_t          n_random_;
      const vector<double>& y_;
   // ----------------------------------------------------------------------
   public:
      // constructor
      mixed_derived(
         size_t n_fixed                    ,
         size_t n_random                   ,
         bool   quasi_fixed                ,
         bool   bool_sparsity              ,
         const CppAD::mixed::d_sparse_rcv&    A_rcv,
         const vector<double>& y           ) :
         cppad_mixed(n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv),
         n_fixed_(n_fixed)     ,
         n_random_(n_random)   ,
         y_(y)
      {  assert( n_fixed == 3);
         assert( y_.size() == n_random_ );
      }
   // ----------------------------------------------------------------------
      // implementation of ran_likelihood
      // Note that theta[2] is not used by ran_likelihood
      template <typename Vector>
      Vector template_ran_likelihood(
         const Vector& theta  ,
         const Vector& u      )
      {  typedef typename Vector::value_type scalar;

         assert( theta.size() == n_fixed_ );
         assert( u.size() == y_.size() );
         Vector vec(1);

         // initialize part of log-density that is always smooth
         vec[0] = scalar(0.0);

         // sqrt_2pi = CppAD::sqrt(8.0 * CppAD::atan(1.0) );

         for(size_t i = 0; i < n_random_; i++)
         {  scalar mu     = u[i] + theta[0];
            scalar sigma  = theta[1];
            scalar res    = (y_[i] - mu) / sigma;

            // p(y_i | u, theta)
            vec[0] += log(sigma) + res * res / scalar(2.0);
            // following term does not depend on fixed or random effects
            // vec[0] += log(sqrt_2pi);

            // p(u_i | theta)
            vec[0] += u[i] * u[i] / scalar(2.0);
            // following term does not depend on fixed or random effects
            // vec[0] += log(sqrt_2pi);
         }
         return vec;
      }
      // a1_vector version of ran_likelihood
      virtual a1_vector ran_likelihood(
         const a1_vector& fixed_vec, const a1_vector& random_vec
      )
      {  return template_ran_likelihood( fixed_vec, random_vec ); }
      // implementation of fix_likelihood
      template <typename Vector>
      Vector template_fix_likelihood(
         const Vector& fixed_vec  )
      {  typedef typename Vector::value_type scalar;

         assert( fixed_vec.size() == n_fixed_ );
         Vector vec(1);

         // initialize part of log-density that is smooth
         vec[0] = scalar(0.0);

         // sqrt_2pi = CppAD::sqrt( 8.0 * CppAD::atan(1.0) );
         for(size_t j = 0; j < n_fixed_; j++)
         {  scalar mu     = scalar(4.0);
            scalar sigma  = scalar(1.0);
            scalar res    = (fixed_vec[j] - mu) / sigma;
            if( j == 2 )
               res = fixed_vec[0] + fixed_vec[1] + fixed_vec[2];
            // This is a Gaussian term, so entire density is smooth
            vec[0]  += res * res / scalar(2.0);
            // following term does not depend on fixed effects
            // vec[0]  += log(sqrt_2pi * sigma);
         }
         return vec;
      }
      // a1_vector version of fix_likelihood
      virtual a1_vector fix_likelihood(const a1_vector& fixed_vec)
      {  return template_fix_likelihood( fixed_vec ); }
   // ----------------------------------------------------------------------
   public:
      // User defined virtual functions
      //
      // ------------------------------------------------------------------
   };
}

bool sample_fixed_2(void)
{
   bool   ok = true;
   double inf = std::numeric_limits<double>::infinity();
   //
   // initialize gsl random number generator
   size_t random_seed = CppAD::mixed::new_gsl_rng(0);
   //
   size_t n_data   = 10;
   size_t n_fixed  = 3;
   size_t n_random = n_data;
   vector<double>
      fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
   fixed_lower[0] = - inf; fixed_in[0] = 2.0; fixed_upper[0] = inf;
   fixed_lower[1] = .01;   fixed_in[1] = 0.5; fixed_upper[1] = inf;
   fixed_lower[2] = - inf; fixed_in[2] = 2.0; fixed_upper[2] = inf;
   //
   // explicit constriants (in addition to l1 terms)
   vector<double> fix_constraint_lower(0), fix_constraint_upper(0);
   //
   vector<double> data(n_data), random_in(n_random);
   for(size_t i = 0; i < n_data; i++)
   {  data[i]       = double(i + 1);
      random_in[i] = 0.0;
   }

   // object that is derived from cppad_mixed
   bool quasi_fixed = true;
   bool bool_sparsity = false;
   CppAD::mixed::d_sparse_rcv A_rcv; // empty matrix
   mixed_derived mixed_object(
      n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, data
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
   std::string random_ipopt_options =
      "Integer print_level               0\n"
      "String  sb                        yes\n"
      "String  derivative_test           second-order\n"
      "Numeric tol                       1e-8\n"
   ;
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
   // check that none of the constraints are active
   // (Note that the Lagragian w.r.t. theta[2] will be zero because
   // it does not affect the objective).
   ok &= solution.fixed_lag.size() == n_fixed;
   for(size_t i = 0; i < n_fixed; i++)
      ok &= solution.fixed_lag[i] == 0.0;
   ok &= solution.fix_con_lag.size() == 0;
   ok &= solution.ran_con_lag.size() == 0;
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
   // Note all entries in information matrix are non-zero
   size_t K = ( (n_fixed + 1) * n_fixed ) / 2;
   ok &= K == information_rcv.nnz();
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
   typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > matrix;
   //
   // compute sample covariance matrix
   matrix sample_cov = matrix::Zero(n_fixed, n_fixed);
   for(size_t i = 0; i < n_sample; i++)
   {  matrix diff(n_fixed, 1);
      for(size_t j = 0; j < n_fixed; j++)
         diff(j, 0) = sample[ i * n_fixed + j] - solution.fixed_opt[j];
      sample_cov += diff * diff.transpose();
   }
   sample_cov *= 1.0 / double(n_sample);
   //
   matrix info_mat = matrix::Zero(n_fixed, n_fixed);
   for(size_t k = 0; k < K; k++)
   {  size_t i = information_rcv.row()[k];
      size_t j = information_rcv.col()[k];
      info_mat(i, j) = information_rcv.val()[k];
      info_mat(j, i) = information_rcv.val()[k];
   }
   matrix cov_mat = info_mat.inverse();
   //
   for(size_t i = 0; i < n_fixed; i++)
   {  for(size_t j = 0; j < n_fixed; j++)
      {  double value = sample_cov(i, j);
         double check = cov_mat(i, j);
         double scale = std::sqrt( cov_mat(i, i) * cov_mat(j, j) );
         ok &= std::fabs(value - check) / scale < .05;
      }
   }
   //
   if( ! ok )
      std::cout << "\nrandom_seed = " << random_seed << "\n";
   //
   CppAD::mixed::free_gsl_rng();
   return ok;
}
// END C++
