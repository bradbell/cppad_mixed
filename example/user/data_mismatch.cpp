// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
$begin data_mismatch.cpp$$
$spell
   CppAD
   cppad
   hes
   eval
   interp
   xam
$$

$section Random Effects Variance May Cause Data Mismatch$$

$head Model$$
$latex \[
   \B{p}( z | \theta ) \sim \B{N} ( \theta , \sigma_z^2 )
\] $$
$latex \[
   \B{p}( y | \theta ) \sim \B{N} [ \exp( u ) \theta , \sigma_y^2 ]
\] $$
$latex \[
   \B{p}( u | \theta ) \sim \B{N} ( 0 , \sigma_u^2 )
\] $$
The fixed likelihood
$cref/g(theta)/theory/Fixed Likelihood, g(theta)/$$
is
$latex \[
g( \theta ) = \frac{1}{2} \left[
   \log ( 2 \pi \sigma_z^2 ) + ( z - \theta )^2
\right]
\] $$
The random likelihood
$cref/f(theta, u)
   /theory/
   Random Likelihood, f(theta, u)
/$$
is
$latex \[
f(\theta , u ) = \frac{1}{2} \left[
   \log ( 2 \pi \sigma_u^2 ) + u^2 / \sigma_u^2
   +
   \log ( 2 \pi \sigma_y^2 ) + [ y - \exp(u) \theta ]^2 / \sigma_y^2
\right]
\] $$

$head Mismatch$$
In the case where $latex y = z$$, one might expect the solution
$latex \theta = z$$ and $latex u = 0$$ because
all the data residuals are zero, $latex \B{p}(u | \theta )$$
is maximal, and there is no prior for $latex \theta$$.
This example demonstrates that $latex \theta = z$$ and $latex u = 0$$
may not be optimal for the this case.
To be specific it shows that the derivative of
$cref/L(theta)/theory/Objective/Fixed Effects Objective, L(theta)/$$
may be non-zero.

$head Theory$$
See the $tref theory$$ section for the
theory behind the calculations below:

$head Derivatives$$
$latex \[
\begin{array}{rcl}
g_\theta ( \theta )
& = &
( \theta - z ) \sigma_z^{-2}
\\
f_\theta ( \theta , u )
& = &
[ \exp(u) \theta - y ] \exp(u) \sigma_y^{-2}
\\
f_u ( \theta , u )
& = &
u / \sigma_u^2 + [ \exp(u) \theta - y ] \exp(u) \theta \sigma_y^{-2}
\\
f_{u,u} ( \theta , u )
& = &
\sigma_u^{-2}
+
[ 2 \exp(u) \theta - y ] \exp(u) \theta \sigma_y^{-2}
\\
f_{u,\theta} ( \theta , u )
& = &
[ 2 \exp(u) \theta - y ] \exp(u) \sigma_y^{-2}
\\
\hat{u}_\theta ( \theta )
& = &
- f_{u,\theta} [ \theta , \hat{u} ( \theta ) ]
/
f_{u,u} [ \theta , \hat{u} ( \theta ) ]
\\
f_{u,u,u} ( \theta , u )
& = &
[ 4 \exp(u) \theta - y ] \exp(u) \theta \sigma_y^{-2}
\\
f_{u,u,\theta} ( \theta , u )
& = &
[ 4 \exp(u) \theta - y ] \exp(u) \sigma_y^{-2}
\end{array}
\] $$

$head Objective$$
$latex \[
\begin{array}{rcl}
h( \theta , u )
& = &
\frac{1}{2} \log f_{u,u} ( \theta , u )
+
f( \theta , u )
-
\log( 2 \pi )
\\
h_\theta ( \theta , u )
& = &
\frac{1}{2} f_{u,u,\theta} ( \theta , u ) / f_{u,u} ( \theta , u )
+
f_\theta ( \theta , u )
\\
h_u ( \theta , u )
& = &
\frac{1}{2} f_{u,u,u} ( \theta , u ) / f_{u,u} ( \theta , u )
+
f_u ( \theta , u )
\\
L( \theta )
& = &
h [ \theta , \hat{u} ( \theta ) ] + g ( \theta )
\\
L_\theta ( \theta )
& = &
h_\theta [ \theta , \hat{u} ( \theta ) ]
+
h_u [ \theta , \hat{u} ( \theta ) ] \hat{u}_\theta ( \theta )
+
g_\theta ( \theta )
\end{array}
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
   using CppAD::log;
   using CppAD::exp;
   using CppAD::AD;
   //
   using CppAD::mixed::d_sparse_rcv;
   using CppAD::mixed::a1_double;
   using CppAD::mixed::d_vector;
   using CppAD::mixed::a1_vector;

   class mixed_derived : public cppad_mixed {
   private:
      const double y_, z_;
      const double sigma_u_, sigma_y_, sigma_z_;
   public:
      // constructor
      mixed_derived(
         size_t                 n_fixed       ,
         size_t                 n_random      ,
         double                 y             ,
         double                 z             ,
         double                 sigma_u       ,
         double                 sigma_y       ,
         double                 sigma_z       ) :
         cppad_mixed(n_fixed, n_random) ,
         y_(y)                          ,
         z_(z)                          ,
         sigma_u_(sigma_u)              ,
         sigma_y_(sigma_y)              ,
         sigma_z_(sigma_z )
      {  assert( n_fixed == 1 );
         assert( n_random == 1 );
      }
      // implementation of fix_likelihood as p(z|theta) * p(theta)
      a1_vector fix_likelihood(
         const a1_vector&         fixed_vec  ) override
      {
         a1_double theta = fixed_vec[0];

         // initialize log-density
         a1_vector vec(1);
         vec[0] = 0.0;

         // compute this factor once
         // sqrt_2pi = CppAD::sqrt( 8.0 * CppAD::atan(1.0) );

         // Data term p(z|theta)
         a1_double res  = (z_ - theta) / sigma_z_;
         vec[0]    += res * res / 2.0;
         // the following term does not depend on fixed effects
         // vec[0]    += log(sqrt_2pi * sigma_z_ );

         // prior term p(theta)

         return vec;
      }
      // ------------------------------------------------------------------
      // implementation of ran_likelihood as p(y|theta, u) * p(u|theta)
      a1_vector ran_likelihood(
         const a1_vector&         fixed_vec  ,
         const a1_vector&         random_vec ) override
      {
         a1_double theta = fixed_vec[0];
         a1_double u     = random_vec[0];

         // initialize log-density
         a1_vector vec(1);
         vec[0] = 0.0;

         // compute this factors once
         // sqrt_2pi = CppAD::sqrt( 8.0 * CppAD::atan(1.0) );

         // Data term p(y|theta,u)
         a1_double res  = (y_ - exp(u) * theta) / sigma_y_;
         vec[0]    += res * res / 2.0;
         // the following term does not depend on fixed or random effects
         // vec[0]    += log(sqrt_2pi * sigma_y_ );

         // prior term p(u|theta)
         res        = u / sigma_u_;
         vec[0]    += res * res / 2.0;
         // the following term does not depend on fixed or random effects
         // vec[0]    += log(sqrt_2pi * sigma_u_ );

         return vec;
      }
      // ==================================================================
      // Routines used to check that objective derivative is zero at solution
      // g(thata)    = (theta - z)^2           / (2.0 * sigma_z * sigma_z)
      // f(theta, u) = ( exp(u) * theta - y)^2 / (2.0 * sigma_y * sigma_y )
      //             + u^2  / (2.0 * sigma_u * sigma_u);
      double g_theta(double theta)
      {  return (theta - z_) / (sigma_z_ * sigma_z_); }
      //
      double f_theta(double theta, double u)
      {  double num_theta    = 2.0 * (exp(u) * theta - y_) * exp(u);
         double ratio_theta  = num_theta / (2.0 * sigma_y_ * sigma_y_);
         return ratio_theta;
      }
      double f_u(double theta, double u)
      {  double num_u    = 2.0 * (exp(u) * theta - y_) * exp(u) * theta;
         double ratio_u  = num_u / (2.0 * sigma_y_ * sigma_y_);
         return ratio_u + u / (sigma_u_ * sigma_u_);
      }
      double f_uu(double theta, double u)
      {  double term     = exp(u) * theta;
         double num_uu   = 2.0 * (term * term + (term - y_) * term);
         double ratio_uu = num_uu / (2.0 * sigma_y_ * sigma_y_);
         return ratio_uu + 1.0 / (sigma_u_ * sigma_u_);
      }
      double f_utheta(double theta, double u)
      {  double num_utheta    = 2.0 * (exp(u) * theta - y_) * exp(u);
         num_utheta          += 2.0 * exp(u) * exp(u) * theta;
         double ratio_utheta  = num_utheta / (2.0 * sigma_y_ * sigma_y_);
         return ratio_utheta;
      }
      double uhat_theta(double theta, double uhat)
      {  double ret = - f_utheta(theta, uhat) / f_uu(theta, uhat);
         return ret;
      }
      double f_uuu(double theta, double u)
      {  double term = exp(u) * theta;
         // term_u = term
         // num_uu = 2.0 * term * (2.0 * term  - y)
         //        = 4.0 * term^2 - 2.0 * term * y
         double num_uuu  = 8.0 * term * term  - 2.0 * term * y_;
         double ratio_uuu = num_uuu / (2.0 * sigma_y_ * sigma_y_);
         return ratio_uuu;
      }
      double f_uutheta(double theta, double u)
      {  double term = exp(u) * theta;
         // term_theta = exp(u)
         // num_uu = 2.0 * term * (2.0 * term  - y)
         //        = 4.0 * term^2 - 2.0 * term * y
         double num_uutheta = 8.0 * term * exp(u) - 2.0 * y_ * exp(u);
         double ratio_uutheta = num_uutheta / (2.0 * sigma_y_ * sigma_y_);
         return ratio_uutheta;
      }
      double h_theta(double theta, double u)
      {  double ret = 0.5 * f_uutheta(theta, u) / f_uu(theta, u);
         ret       += f_theta(theta, u);
         return ret;
      }
      double h_u(double theta, double u)
      {  double ret = 0.5 * f_uuu(theta, u) / f_uu(theta, u);
         ret       += f_u(theta, u);
         return ret;
      }
      double L_theta(double theta , double uhat)
      {  double ret = h_theta(theta, uhat);
         ret       += h_u(theta, uhat) * uhat_theta(theta, uhat);
         ret       += g_theta(theta);
         return ret;
      }
   };

}

bool data_mismatch_xam(void)
{
   bool   ok = true;
   double inf = std::numeric_limits<double>::infinity();

   size_t n_fixed  = 1;
   size_t n_random = 1;
   double z        = 0.05;
   double y        = z;
   double sigma_u  = 0.1;
   double sigma_y  = 0.1 * y;
   double sigma_z  = 0.1 * z;
   d_vector fixed_in(n_fixed), random_in(n_random);
   d_vector fixed_lower(n_fixed), fixed_upper(n_fixed);
   fixed_lower[0] = -inf;
   fixed_in[0]    = z;
   fixed_upper[0] = + inf;
   random_in[0]   = 0.0;
   //
   // no constriants
   d_vector fix_constraint_lower(0), fix_constraint_upper(0);
   //
   //
   // object that is derived from cppad_mixed
   mixed_derived mixed_object(
      n_fixed, n_random,
      y, z, sigma_u, sigma_y, sigma_z

   );
   mixed_object.initialize(fixed_in, random_in);
   //
   // compute the derivative of the objective at the starting point
   double theta_in   = fixed_in[0];
   double u_in       = random_in[0];
   double L_theta_in = mixed_object.L_theta(theta_in, u_in);
   //
   // derivative corresponding to theta = z = y is greater than 1
   ok &= fabs( L_theta_in ) > 1.0;
   //
   // optimize the fixed effects using quasi-Newton method
   std::string fixed_ipopt_options =
      "Integer print_level               0\n"
      "String  sb                        yes\n"
      "String  derivative_test           adaptive\n"
      "Numeric tol                       1e-8\n"
      "Integer max_iter                  15\n"
   ;
   std::string random_ipopt_options =
      "Integer print_level               0\n"
      "String  sb                        yes\n"
      "String  derivative_test           second-order\n"
      "Numeric tol                       1e-8\n"
   ;
   //
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
   d_vector random_out = mixed_object.optimize_random(
      random_ipopt_options, fixed_out, random_lower, random_upper, random_in
   );
   //
   // compute the derivative of the objective at the final point
   double theta_out   = fixed_out[0];
   double u_out       = random_out[0];
   double L_theta_out = mixed_object.L_theta(theta_out, u_out);
   //
   // derivative corresponding to solution is less that 1e-7
   ok &= fabs( L_theta_out ) <= 1e-7;
   //
   // Now demonstrate that the solution is still close to the expected values
   ok &= fabs( theta_out / z - 1.0 ) <= 1e-2;
   ok &= fabs( u_out ) <= 1e-2;
   //
   return ok;
}
// END C++
