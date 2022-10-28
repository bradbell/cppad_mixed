// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
// case with no cross terms; i.e, f_{u,theta} ( theta , u ) is zero

# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>

namespace {
   using CppAD::vector;
   using CppAD::AD;
   using CppAD::mixed::d_sparse_rcv;
   //
   using CppAD::mixed::a1_double;

   class mixed_derived : public cppad_mixed {
   private:
      const vector<double>& y_;
   public:
      // constructor
      mixed_derived(
         size_t n_fixed                    ,
         size_t n_random                   ,
         const CppAD::mixed::d_sparse_rcv&    A_rcv,
         const vector<double>& y           )
         :
         // quasi_fixed = false, bool_sparsity = true
         cppad_mixed(n_fixed, n_random, false, true, A_rcv) ,
         y_(y)
      {  assert( n_fixed == 2);
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

         for(size_t i = 0; i < y_.size(); i++)
         {  scalar mu     = u[i];
            scalar sigma  = scalar(0.5);
            scalar res    = (y_[i] - mu) / sigma;

            // p(y_i | u, theta)
            vec[0] += CppAD::log(sqrt_2pi * sigma) + res*res / scalar(2.0);

            // p(u_i | theta)
            vec[0] += CppAD::log(sqrt_2pi) + u[i] * u[i] / scalar(2.0);
         }
         return vec;
      }
      // a1_vector version of ran_likelihood
      virtual a1_vector ran_likelihood(
         const a1_vector& fixed_vec, const a1_vector& random_vec
      )
      {  return template_ran_likelihood( fixed_vec, random_vec ); }
   public:
      //
   };
}

bool ran_obj_tst(void)
{
   bool   ok = true;
   double eps = 100. * std::numeric_limits<double>::epsilon();

   size_t n_data   = 1;
   size_t n_fixed  = 2;
   size_t n_random = n_data;
   vector<double> data(n_data), fixed_vec(n_fixed), random_vec(n_random);
   vector<double> uhat(n_random);

   fixed_vec[0] = 2.0;
   fixed_vec[1] = 0.5;
   for(size_t i = 0; i < n_data; i++)
   {  data[i]       = double(i + 1);
      random_vec[i] = double(i) / double(n_data);
   }

   // object that is derived from cppad_mixed
   CppAD::mixed::d_sparse_rcv A_rcv; // empty matrix
   mixed_derived mixed_object(n_fixed, n_random, A_rcv, data);
   mixed_object.initialize(fixed_vec, random_vec);

   // lower and upper limits for random effects
   double inf = std::numeric_limits<double>::infinity();
   vector<double> random_lower(n_random), random_upper(n_random);
   for(size_t i = 0; i < n_random; i++)
   {  random_lower[i] = -inf;
      random_upper[i] = +inf;
   }

   // optimize the random effects
   std::string options;
   options += "Integer print_level 0\n";
   options += "String  sb          yes\n";
   options += "String  derivative_test second-order\n";
   uhat = mixed_object.optimize_random(
      options, fixed_vec, random_lower, random_upper, random_vec
   );

   // factor f_{u,u} ( theta , uhat )
   mixed_object.update_factor(fixed_vec, uhat);

   // compute total derivative of random part of objective
   // vector<double> r_fixed(n_fixed);
   // mixed_object.ran_obj_jac(fixed_vec, uhat, r_fixed);


   // For this case the Laplace approximation is exactly equal the integral
   // p(y_i | theta ) = integral of p(y_i | theta , u) p(u | theta) du
   // Furthermore p(y_i | theta ) is N( theta[0], 1 + 0.5^2 )

   // check the random part of the objective
   double r        = mixed_object.ran_obj_eval(fixed_vec, uhat);
   double check    = 0.0;
   double sqrt_2pi = CppAD::sqrt(8.0 * CppAD::atan(1.0) );
   double sigma    = CppAD::sqrt( 1.0 + 0.5 * 0.5 );
   for(size_t i = 0; i < n_data; i++)
   {  double res   = data[i] / sigma;
      check       += CppAD::log(sqrt_2pi * sigma);
      check       += res * res / 2.0;
   }
   ok &= fabs( r / check - 1.0 ) <= eps;

   // check jacobian of random part of objective
   vector<double> r_fixed(n_fixed);
   mixed_object.ran_obj_jac(fixed_vec, uhat, r_fixed);
   //
   ok &= fabs( r_fixed[0] - 0.0 ) <= eps;
   ok &= fabs( r_fixed[1] - 0.0 ) <= eps;

   return ok;
}
