// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
// information matrix and sample_fixed in case where with no random effects.

# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/sparse_mat_info.hpp>
# include <cppad/mixed/manage_gsl_rng.hpp>

namespace {
   using CppAD::vector;
   using CppAD::log;
   using CppAD::AD;
   using CppAD::mixed::d_sparse_rcv;

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
      {  assert( y_.size() == n_fixed_ );
         assert( n_random_ == 0 );
      }
   // ----------------------------------------------------------------------
   private:
      template <class Float>
      vector<Float> implement_fix_likelihood(
         const vector<Float>& theta  )
      {  assert( theta.size() == n_fixed_ );
         assert( theta.size() == y_.size() );
         vector<Float> vec(1);

         // initialize part of log-density that is always smooth
         vec[0] = Float(0.0);

         // pi
         Float sqrt_2pi = Float( CppAD::sqrt(8.0 * CppAD::atan(1.0) ) );

         for(size_t i = 0; i < n_fixed_; i++)
         {  Float mu     = theta[i];
            Float sigma  = double(i + 1);
            Float res    = (y_[i] - mu) / sigma;

            // p(y_i | theta)
            vec[0] += log(sqrt_2pi * sigma) + res * res / Float(2.0);

         }
         return vec;
      }
   // ----------------------------------------------------------------------
   public:
      virtual vector<a1_double> fix_likelihood(
         const vector<a1_double>& fixed_vec  )
      {  return implement_fix_likelihood(fixed_vec); }
      // ------------------------------------------------------------------
   };
}

bool no_random_info(void)
{
   bool   ok = true;
   double inf = std::numeric_limits<double>::infinity();
   double eps = 10. * std::numeric_limits<double>::epsilon();
   //
   // initialize gsl random number generator
   size_t random_seed = CppAD::mixed::new_gsl_rng(0);
   //
   size_t n_data   = 10;
   size_t n_fixed  = n_data;
   size_t n_random = 0;
   vector<double>
      fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
   for(size_t i = 0; i < n_fixed; i++)
   {  fixed_lower[i] = -inf;
      fixed_in[i]    = 0.0;
      fixed_upper[i] = + inf;
   }
   //
   // set one of the fixed effect to be zero
   size_t set_zero_index = 2;
   fixed_lower[set_zero_index] = 0.0;
   fixed_in[set_zero_index]    = 0.0;
   fixed_upper[set_zero_index] = 0.0;
   //
   // explicit constriants (in addition to l1 terms)
   vector<double> fix_constraint_lower(0), fix_constraint_upper(0);
   //
   vector<double> data(n_data);
   for(size_t i = 0; i < n_data; i++)
      data[i]       = double(i + 1);

   // zero length vecors
   vector<double> random_lower, random_in, random_upper, random_opt;

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
      "String  derivative_test           first-order\n"
      "String  derivative_test_print_all yes\n"
      "Numeric tol                       1e-8\n"
      "Integer max_iter                  50\n"
   ;
   std::string random_ipopt_options =
      "Integer print_level     0\n"
      "String  sb              yes\n"
      "String  derivative_test second-order\n"
      "Numeric tol             1e-8\n"
   ;
   //
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
   // compute corresponding information matrix
   d_sparse_rcv
   information_rcv = mixed_object.information_mat(solution, random_opt);
   //
   // The infromation matrix is diagonal
   ok  &= information_rcv.nnz() == n_fixed;
   //
   // check Hessian values
   const CppAD::mixed::s_vector& row( information_rcv.row() );
   const CppAD::mixed::s_vector& col( information_rcv.col() );
   const CppAD::mixed::d_vector& val( information_rcv.val() );
   for(size_t i = 0; i < n_fixed; i++)
   {  ok &= row[i] == i;
      ok &= col[i] == i;
      double check = 1.0 / double( (i + 1) * (i + 1) );
      ok &= CppAD::NearEqual( val[i], check, eps, eps);
   }
   //
   // make sure can sample from fixed effects
   size_t sample_size = 1;
   vector<double> sample(n_fixed * sample_size);
   mixed_object.sample_fixed(
      sample,
      information_rcv,
      solution,
      fixed_lower,
      fixed_upper
   );
   //
   // check that the sample is reasonable
   for(size_t i = 0; i < n_fixed; i++)
   {  double sigma = double(i + 1);
      // note hessian is 1 / sigma^2 along diagonal
      ok   &= std::fabs( sample[i] - solution.fixed_opt[i] ) < 4 * sigma;
      if( i == set_zero_index )
      {  ok &= solution.fixed_opt[i] == 0.0;
         ok &= sample[i] == 0.0;
      }
   }
   if( ! ok )
      std::cout << "\nrandom_seed = " << random_seed << "\n";
   //
   CppAD::mixed::free_gsl_rng();
   return ok;
}
// END C++
