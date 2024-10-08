// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-24 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
@begin zero_random_two.cpp@@
$spell
   CppAD
   cppad
   hes
   eval
   interp
   xam
$$

$section Case Where Optimum is Non-Zero Fixed and Zero Random Effects$$

$head Model$$
$latex \[
   \B{p}( y_i | \theta , u ) \sim \B{N} left( \theta_0 , \theta_1^2 \right)
\] $$
$latex \[
   \B{p}( u_i | \theta ) \sim \B{N} ( 0 , 1 )
\] $$
$latex \[
   \B{p}( \theta ) \sim \B{N} ( 1 , 1 )
\] $$
It follows that the Laplace approximation is exact and
$latex \[
   \B{p}( y_i | \theta ) \sim \B{N} \left( \theta_0 , \theta_1^2 \right)
\] $$
The corresponding objective for the fixed effects is equivalent to:
$latex \[
F( \theta ) = \frac{1}{2} \left[
   ( \theta_0 - 1 )^2 + ( \theta_1 - 1 )^2 +
      N \log \left( \theta_1^2 \right) +
      \theta_1^{-2} \sum_{i=0}^{N-1} ( y_i - \theta_0 )^2
\right]
\] $$
The constraints on the fixed effect are
$latex \[
   - \infty \leq \theta_0 \leq + \infty
   \R{\; and \;}
   0.1 \leq \theta_1 \leq 100
\] $$

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
   using CppAD::mixed::d_sparse_rcv;
   //
   using CppAD::mixed::a1_double;

   class mixed_derived : public cppad_mixed {
   private:
      const size_t n_fixed_;
# ifndef NDEBUG
      const size_t n_random_;
# endif
      const vector<double> y_;
   // ----------------------------------------------------------------------
   public:
      // constructor
      mixed_derived(
         size_t n_fixed                    ,
         size_t n_random                   ,
         const CppAD::mixed::d_sparse_rcv&    A_rcv,
         const vector<double>& y           ) :
         // quasi_fixed = false, bool_sparsity = true
         cppad_mixed(n_fixed, n_random, false, true, A_rcv) ,
         n_fixed_(n_fixed)                          ,
# ifndef NDEBUG
         n_random_(n_random)                        ,
# endif
         y_(y)
      {  assert( n_fixed_ == 2);}
   // ----------------------------------------------------------------------
      // implementation of ran_likelihood
      template <typename Vector>
      Vector template_ran_likelihood(
         const Vector& theta  ,
         const Vector& u      )
      {  typedef typename Vector::value_type scalar;

         assert( theta.size() == n_fixed_ );
         assert( u.size()     == n_random_ );
         Vector vec(1);

         // initialize part of log-density that is always smooth
         vec[0] = scalar(0.0);

         // square root of 2 * pi
         scalar sqrt_2pi = scalar(
             CppAD::sqrt(8.0 * CppAD::atan(1.0)
         ));

         for(size_t i = 0; i < u.size(); i++)
         {  // p(u_i | theta)
            vec[0] += log(sqrt_2pi) + u[i] * u[i] / scalar(2.0);
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

         Vector vec(1);

         // initialize part of log-density that is smooth
         vec[0] = scalar(0.0);

         // compute these factors once
         scalar sqrt_2pi = scalar(
             CppAD::sqrt( 8.0 * CppAD::atan(1.0)
         ));

         // prior for fixed effects
         for(size_t j = 0; j < n_fixed_; j++)
         {  scalar mu     = scalar(1.0);
            scalar sigma  = scalar(1.0);
            scalar res    = (fixed_vec[j] - mu) / sigma;

            // p(theta_j )
            vec[0]  += log(sqrt_2pi * sigma) + res * res / scalar(2.0);
         }


         // data given fixed effects
         for(size_t i = 0; i < y_.size(); i++)
         {  scalar mu     = fixed_vec[0];
            scalar sigma  = fixed_vec[1];
            scalar res    = (y_[i] - mu) / sigma;

            // p(y_i | theta)
            vec[0] += log(sqrt_2pi * sigma) + res * res / scalar(2.0);
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
      dF[0]  = (theta[0] - 1.0) + sum / (theta[1] * theta[1]);
      dF[1]  = theta[1] - 1.0;
      dF[1] += double(data.size()) / theta[1];
      dF[1] -= sumsq  / (theta[1] * theta[1] * theta[1]);
      //
      return dF;
   }
}

bool zero_random_two(void)
{
   bool   ok = true;
   double inf = std::numeric_limits<double>::infinity();
   double tol = 1e-8;

   size_t n_data   = 10;
   size_t n_fixed  = 2;
   size_t n_random = 3;;
   vector<double>
      fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
   fixed_lower[0] = - inf; fixed_in[0] = 2.0; fixed_upper[0] = inf;
   fixed_lower[1] = .01;   fixed_in[1] = 0.5; fixed_upper[1] = inf;
   //
   // explicit constriants (in addition to l1 terms)
   vector<double> fix_constraint_lower(0), fix_constraint_upper(0);
   //
   vector<double> data(n_data);
   for(size_t i = 0; i < n_data; i++)
      data[i]       = double(i + 1);
   //
   vector<double>  random_in(n_random);
   for(size_t i = 0; i < n_random; i++)
      random_in[i] = 1.0;

   // object that is derived from cppad_mixed
   CppAD::mixed::d_sparse_rcv A_rcv; // empty matrix
   mixed_derived mixed_object(n_fixed, n_random, A_rcv, data);
   mixed_object.initialize(fixed_in, random_in);

   // optimize the fixed effects
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

   // Compute the optimal random effects
   vector<double> random_out = mixed_object.optimize_random(
      random_ipopt_options, fixed_out, random_lower, random_upper, random_in
   );
   // partial of p(u | theta) w.r.t u_i is equal to u_i
   for(size_t i = 0; i < n_random; i++)
      ok &= fabs( random_out[i] ) <= 5.0 * tol;

   return ok;
}
// END C++
