// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-24 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
$section Optimize Fixed Effects With No fix_likelihood: Example and Test$$

$head Model$$
$latex \[
   \B{p}( y_i | \theta , u ) \sim \B{N} ( u_i + \theta_0 , \theta_1^2 )
\] $$
$latex \[
   \B{p}( u_i | \theta ) \sim \B{N} ( 0 , 1 )
\] $$
It follows that the Laplace approximation is exact and
$latex \[
   \B{p}( y_i | \theta ) \sim \B{N} \left( \theta_0 , 1 + \theta_1^2 \right)
\] $$
The corresponding objective for the fixed effects is equivalent to:
$latex \[
F( \theta ) = \frac{1}{2} \left[
      N \log \left( 1 + \theta_1^2 \right) +
      ( 1 + \theta_1^2)^{-1} \sum_{i=0}^{N-1} ( y_i - \theta_0 )^2
\right]
\] $$
The constraints on the fixed effect are
$latex \[
   - \infty \leq \theta_0 \leq + \infty
   \R{\; and \;}
   0.1 \leq \theta_1 \leq 100
\] $$

*/
// BEGIN C++
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
# ifndef NDEBUG
      const size_t          n_fixed_;
# endif
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
# ifndef NDEBUG
         n_fixed_(n_fixed)     ,
# endif
         n_random_(n_random)   ,
         y_(y)
      {  assert( n_fixed == 2);
         assert( y_.size() == n_random_ );
      }
   // ----------------------------------------------------------------------
      // implementation of ran_likelihood
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

         // pi
         scalar sqrt_2pi = scalar(
             CppAD::sqrt(8.0 * CppAD::atan(1.0)
         ));

         for(size_t i = 0; i < n_random_; i++)
         {  scalar mu     = u[i] + theta[0];
            scalar sigma  = theta[1];
            scalar res    = (y_[i] - mu) / sigma;

            // p(y_i | u, theta)
            vec[0] += log(sqrt_2pi * sigma) + res * res / scalar(2.0);

            // p(u_i | theta)
            vec[0] += log(sqrt_2pi) + u[i] * u[i] / scalar(2.0);
         }
         return vec;
      }
      // a1_vector version of ran_likelihood
      virtual a1_vector ran_likelihood(
         const a1_vector& fixed_vec, const a1_vector& random_vec
      )
      {  return template_ran_likelihood( fixed_vec, random_vec ); }
   // ----------------------------------------------------------------------
   public:
      // User defined virtual functions
      //
   };
   // derivative of objective
   vector<double> objective_fixed(
      const vector<double>& data   ,
      const vector<double>& theta  )
   {  vector<double> dF(2);
      //
      // compute partials of F
      double sum   = 0.0;
      double sumsq = 0.0;
      for(size_t i = 0; i < data.size(); i++)
      {  sum   += theta[0] - data[i];
         sumsq += (theta[0] - data[i]) * (theta[0] - data[i]);
      }
      double den = 1.0 + theta[1] * theta[1];
      dF[0]  = sum / den;
      dF[1]  = 0.0;
      dF[1] += double(data.size()) * theta[1] / den;
      dF[1] -= sumsq * theta[1]  / (den * den);
      //
      return dF;
   }
}

bool no_fix_likelihood(void)
{
   bool   ok = true;
   double inf = std::numeric_limits<double>::infinity();
   double tol = 1e-8;

   size_t n_data   = 10;
   size_t n_fixed  = 2;
   size_t n_random = n_data;
   vector<double>
      fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
   fixed_lower[0] = - inf; fixed_in[0] = 2.0; fixed_upper[0] = inf;
   fixed_lower[1] = .01;   fixed_in[1] = 0.5; fixed_upper[1] = inf;
   //
   // explicit constriants (in addition to l1 terms)
   vector<double> fix_constraint_lower(0), fix_constraint_upper(0);
   //
   vector<double> data(n_data), random_in(n_random);
   for(size_t i = 0; i < n_data; i++)
   {  data[i]       = double(i + 1);
      random_in[i] = 0.0;
   }

   // use full Newton method to test Hessian calcualtion
   // when there is no fix_likelihood
   bool quasi_fixed   = false;
   bool bool_sparsity = true;

   // object that is derived from cppad_mixed
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
      "Integer print_level     0\n"
      "String  sb              yes\n"
      "String  derivative_test second-order\n"
      "Numeric tol             1e-8\n"
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

   // deriative of objective at fixed_in and fixed_out
   vector<double> dF_in  = objective_fixed(data, fixed_in);
   vector<double> dF_out = objective_fixed(data, fixed_out);

   // scaling for objective
   double scale = std::max( std::fabs( dF_in[0] ), std::fabs( dF_in[1] ) );
   scale = 1.0 / scale;

   // Note that no constraints are active, (not even the l1 terms)
   // so the partials should be zero.
   ok &= fabs( scale * dF_out[0] ) <= tol;
   ok &= fabs( scale * dF_out[1] ) <= tol;

   return ok;
}
// END C++
