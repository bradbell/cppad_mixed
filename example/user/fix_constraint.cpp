// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-24 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin fix_constraint.cpp}

Using Constraints: Example and Test
###################################

Model
*****

.. math::

   \B{p}( y_i | \theta , u ) \sim \B{N} ( u_i + \theta_i , 1 )

.. math::

   \B{p}( u_i | \theta ) \sim \B{N} ( 0 , 1 )

.. math::

   \B{p}( \theta ) \sim \B{U} ( - \infty , + \infty )

where :math:`\B{U} ( - \infty ,  + \infty )` is the improper uniform prior
on :math:`[- \infty , + \infty ]`.
It follows that the Laplace approximation is exact and

.. math::

   \B{p}( y_i | \theta ) \sim \B{N} ( \theta_i , 2 )

The corresponding objective for the fixed effects is equivalent to:

.. math::

   \frac{1}{2} \sum_{i=0}^{N-1} ( y_i - \theta_i )^2

For this problem we add the explicit constraint

.. math::

   \frac{1}{2} \sum_i \theta_i^2 = 1;

The corresponding Lagrangian is

.. math::

   L( \theta , \lambda ) =
      \frac{1}{2} \sum_{i=0}^{N-1} ( y_i - \theta_i )^2
         + \lambda \left( \frac{1}{2} \sum_i \theta_i^2 - 1 \right)

{xrst_literal
   // BEGIN C++
   // END C++
}

{xrst_end fix_constraint.cpp}
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
      const d_vector&       y_;
   public:
      // constructor
      mixed_derived(
         size_t                   n_fixed        ,
         size_t                   n_random       ,
         const d_vector&          y              ) :
         cppad_mixed(n_fixed, n_random) ,
         y_(y)
      {}
      // implementation of ran_likelihood
      a1_vector ran_likelihood(
         const a1_vector&         theta  ,
         const a1_vector&         u      ) override
      {
         assert( u.size() == y_.size() );
         assert( theta.size() == y_.size() );
         a1_vector vec(1);

         // initialize part of log-density that is always smooth
         vec[0] = 0.0;

         // pi
         // sqrt_2pi = CppAD::sqrt(8.0 * CppAD::atan(1.0) );

         for(size_t i = 0; i < y_.size(); i++)
         {  a1_double mu     = u[i] + theta[i];
            a1_double sigma  = 1.0;
            a1_double res    = (y_[i] - mu) / sigma;

            // p(y_i | u, theta)
            vec[0] += res*res / 2.0;
            // following term does not depend on fixed or random effects
            // vec[0] += log(sqrt_2pi * sigma);

            // p(u_i | theta)
            vec[0] += u[i] * u[i] / 2.0;
            // following term does not depend on fixed or random effects
            // vec[0] += log(sqrt_2pi);
         }
         return vec;
      }
      // ------------------------------------------------------------------
      // fix_constraint
      a1_vector template_fix_constraint(
         const a1_vector&         fixed_vec  )
      {
         a1_vector ret_val(1);
         //
         ret_val[0] = 0.0;
         for(size_t i = 0; i < fixed_vec.size(); i++)
            ret_val[0] += fixed_vec[i] * fixed_vec[i];
         ret_val[0] /= 2.0;
         //
         return ret_val;
      }
      // a1_vector version of fix_constraint
      a1_vector fix_constraint(const a1_vector& fixed_vec) override
      {  return template_fix_constraint( fixed_vec ); }
   };
}

bool fix_constraint_xam(void)
{
   bool   ok = true;
   double inf = std::numeric_limits<double>::infinity();
   double tol = 1e-8;

   size_t n_data   = 3;
   size_t n_fixed  = n_data;
   size_t n_random = n_data;
   d_vector
      fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
   for(size_t i = 0; i < n_fixed; i++)
   {  fixed_lower[i] = - inf;
      fixed_in[i]    = 0.1;
      fixed_upper[i] = inf;
   }
   //
   // explicit constriants (in addition to l1 terms)
   d_vector fix_constraint_lower(1), fix_constraint_upper(1);
   fix_constraint_lower[0] = 1.0;
   fix_constraint_upper[0] = 1.0;
   //
   d_vector data(n_data), random_in(n_random);
   for(size_t i = 0; i < n_data; i++)
   {  data[i]       = double(i + 1);
      random_in[i] = 0.0;
   }

   // object that is derived from cppad_mixed
   mixed_derived mixed_object(n_fixed, n_random, data);
   mixed_object.initialize(fixed_in, random_in);

   // optimize the fixed effects using full Newton method
   std::string fixed_ipopt_options =
      "Integer print_level               0\n"
      "String  sb                        yes\n"
      "String  derivative_test           adaptive\n"
      "String  derivative_test_print_all yes\n"
      "Numeric tol                       1e-8\n"
      "Integer max_iter                  15\n"
   ;
   // random_ipopt_options is non-empty, so using ipopt for random effects
   std::string random_ipopt_options =
      "Integer print_level 0\n"
      "String  sb          yes\n"
      "String  derivative_test second-order\n"
   ;
   // lower and upper limits for random effects
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
   // check constraint
   double sum = 0.0;
   for(size_t i = 0; i < n_fixed; i++)
      sum += fixed_out[i] * fixed_out[i];
   ok &= fabs( sum / 2.0 - 1.0 ) <= tol;

   // compute lagranges multiplier by averaging
   sum = 0.0;
   for(size_t i = 0; i < n_fixed; i++)
      sum += (fixed_out[i] - data[i]) / fixed_out[i];
   double lambda = sum / double(n_fixed);

   // check partials of Lagragian w.r.t fixed effects
   for(size_t i = 0; i < n_fixed; i++)
   {  double err  = data[i] - fixed_out[i] + lambda * fixed_out[i];
      ok         &= fabs(err) < tol;
   }
   return ok;
}
// END C++
