// $Id$
/* --------------------------------------------------------------------------
dismod_at: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin der_var_hes.cpp$$
$spell
	cppad
	hes
	eval
	interp
	xam
$$

$section C++ cppad_mixed: Derivative Test where Hessian Depends on Random Effects$$

$head Model$$
$latex \[
	\B{p}( y | \theta ) \sim \B{N} [ \exp( u ) \theta , \sigma_y^2 ]
\] $$
$latex \[
	\B{p}( u | \theta ) \sim \B{N} ( 0 , \sigma_u^2 )
\] $$
The negative log-likelihood for the random effects
$cref/f(theta, u)
	/cppad_mixed_theory/
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

$head Theory$$
See the $tref cppad_mixed_new_theory$$ section for the
theory behind the calculations below:

$head Derivatives$$
$latex \[
\begin{array}{rcl}
f_\theta^{(1)} ( \theta , u )
& = &
[ \exp(u) \theta - y ] \exp(u) \sigma_y^{-2}
\\
f_u^{(1)} ( \theta , u )
& = &
u / \sigma_u^2 + [ \exp(u) \theta - y ] \exp(u) \theta \sigma_y^{-2}
\\
f_{uu}^{(2)} ( \theta , u )
& = &
\sigma_u^{-2}
+
[ 2 \exp(u) \theta - y ] \exp(u) \theta \sigma_y^{-2}
\\
f_{u \theta}^{(2)} ( \theta , u )
& = &
[ 2 \exp(u) \theta - y ] \exp(u) \sigma_y^{-2}
\\
\hat{u}_\theta^{(1)} ( \theta )
& = &
- f_{u \theta}^{(2)} [ \theta , \hat{u} ( \theta ) ]
/
f_{uu}^{(2)} [ \theta , \hat{u} ( \theta ) ]
\\
f_{uuu}^{(3)} ( \theta , u )
& = &
[ 4 \exp(u) \theta - y ] \exp(u) \theta \sigma_y^{-2}
\\
f_{uu \theta}^{(3)} ( \theta , u )
& = &
[ 4 \exp(u) \theta - y ] \exp(u) \sigma_y^{-2}
\end{array}
\] $$

$head Objective$$
$latex \[
\begin{array}{rcl}
h( \theta , u )
& = &
\frac{1}{2} \log f_{uu}^{(2)} ( \theta , u )
+
f( \theta , u )
-
\log( 2 \pi )
\\
h_\theta^{(1)} ( \theta , u )
& = &
\frac{1}{2} f_{uu \theta}^{(3)} ( \theta , u ) / f_{uu}^{(2)} ( \theta , u )
+
f_\theta^{(1)} ( \theta , u )
\\
h_u^{(1)} ( \theta , u )
& = &
\frac{1}{2} f_{uuu}^{(3)} ( \theta , u ) / f_{uu}^{(2)} ( \theta , u )
+
f_u^{(1)} ( \theta , u )
\\
r( \theta )
& = &
h [ \theta , \hat{u} ( \theta ) ]
\\
r_\theta^{(1)} ( \theta )
& = &
h_\theta^{(1)} [ \theta , \hat{u} ( \theta ) ]
+
h_u^{(1)} [ \theta , \hat{u} ( \theta ) ] \hat{u}^{(1)} ( \theta )
\end{array}
\] $$

$code
$verbatim%test/devel/cppad_mixed/der_var_hes.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <dismod_at/cppad_mixed.hpp>

namespace {
	using CppAD::vector;
	using CppAD::log;
	using CppAD::exp;
	using CppAD::abs;
	using CppAD::AD;

	class mixed_derived : public dismod_at::cppad_mixed {
	private:
		const double y_;
		const double sigma_u_, sigma_y_;
	public:
		// constructor
		mixed_derived(
			size_t n_fixed      ,
			size_t n_random     ,
			double y            ,
			double sigma_u      ,
			double sigma_y      ) :
			// quasi_fixed = false
			dismod_at::cppad_mixed(n_fixed, n_random, false) ,
			y_(y) , sigma_u_(sigma_u), sigma_y_(sigma_y)
		{	assert( n_fixed == 1 );
			assert( n_random == 1 );
		}
	public:
		// ------------------------------------------------------------------
		// implementation of ran_like as p(y|theta, u) * p(u|theta)
		template <class Float>
		vector<Float> implement_ran_like(
			const vector<Float>& fixed_vec  ,
			const vector<Float>& random_vec )
		{	Float theta = fixed_vec[0];
			Float u     = random_vec[0];

			// initialize log-density
			vector<Float> vec(1);
			vec[0] = Float(0.0);

			// compute this factors once
			Float sqrt_2pi = Float( CppAD::sqrt( 8.0 * CppAD::atan(1.0) ) );

			// Data term
			Float res  = (y_ - exp(u) * theta) / sigma_y_;
			vec[0]    += log(sqrt_2pi * sigma_y_ ) + res * res / Float(2.0);

			// prior for u
			res        = u / sigma_u_;
			vec[0]    += log(sqrt_2pi * sigma_u_ ) + res * res / Float(2.0);

			return vec;
		}
		//
		virtual vector<a1_double> fix_like(
			const vector<a1_double>& fixed_vec  )
		{	return a1d_vector(0); } // empty vector
		//
		virtual vector<a2_double> ran_like(
			const vector<a2_double>& fixed_vec   ,
			const vector<a2_double>& random_vec  )
		{	return implement_ran_like(fixed_vec, random_vec); }
		virtual vector<a1_double> ran_like(
			const vector<a1_double>& fixed_vec   ,
			const vector<a1_double>& random_vec  )
		{	return implement_ran_like(fixed_vec, random_vec); }
		//
		virtual vector<a1_double> constraint(
			const vector<a1_double>& fixed_vec  )
		{	return a1d_vector(0); } // empty vector
		//
		virtual void fatal_error(const std::string& error_message)
		{	std::cerr << "Error: " << error_message << std::endl;
			assert(false);
		}
		//
		virtual void warning(const std::string& warning_message)
		{	std::cerr << "Warning: " << warning_message << std::endl;
			assert(false);
		}
		// ==================================================================
		// Routines used to check that objective derivative
		double f_theta(double theta, double u)
		{	double ret = (exp(u) * theta - y_) * exp(u);
			ret        = ret / (sigma_y_ * sigma_y_);
			return ret;
		}
		double f_u(double theta, double u)
		{	double term = exp(u) * theta;
			double ret  = u / (sigma_u_ * sigma_u_);
			ret        += (term - y_) * term / (sigma_y_ * sigma_y_);
			return ret;
		}
		double f_uu(double theta, double u)
		{	double term = exp(u) * theta;
			double ret  = 1.0 / (sigma_u_ * sigma_u_);
			ret        += (2.0 * term - y_) * term / (sigma_y_ * sigma_y_);
			return ret;
		}
		double f_utheta(double theta, double u)
		{	double ret = (2.0 * exp(u) * theta - y_) * exp(u);
			ret        = ret / (sigma_y_ * sigma_y_);
			return ret;
		}
		double uhat_theta(double theta, double uhat)
		{	double ret = - f_utheta(theta, uhat) / f_uu(theta, uhat);
			return ret;
		}
		double f_uuu(double theta, double u)
		{	double term = exp(u) * theta;
			double ret  = (4.0 * term - y_) * term / (sigma_y_ * sigma_y_);
			return ret;
		}
		double f_uutheta(double theta, double u)
		{	double ret = (4.0 * exp(u) * theta - y_) * exp(u);
			ret        = ret / (sigma_y_ * sigma_y_);
			return ret;
		}
		double h_theta(double theta, double u)
		{	double ret = 0.5 * f_uutheta(theta, u) / f_uu(theta, u);
			ret       += f_theta(theta, u);
			return ret;
		}
		double h_u(double theta, double u)
		{	double ret = 0.5 * f_uuu(theta, u) / f_uu(theta, u);
			ret       += f_u(theta, u);
			return ret;
		}
		double r_theta(double theta , double uhat)
		{	double ret = h_theta(theta, uhat);
			ret       += h_u(theta, uhat) * uhat_theta(theta, uhat);
			return ret;
		}
	};

}

bool der_var_hes(void)
{
	bool   ok = true;
	double inf = std::numeric_limits<double>::infinity();
	double eps = 100.0 * std::numeric_limits<double>::epsilon();

	size_t n_fixed  = 1;
	size_t n_random = 1;
	double y        = 0.05;
	double sigma_u  = 0.1;
	double sigma_y  = 0.1 * y;
	vector<double> fixed_vec(n_fixed), random_vec(n_random);
	vector<double> fixed_lower(n_fixed), fixed_upper(n_fixed);
	fixed_lower[0] = -inf;
	fixed_vec[0]   = 1.1 * y;
	fixed_upper[0] = + inf;
	random_vec[0]   = 0.0;
	//
	// no constriants
	vector<double> constraint_lower(0), constraint_upper(0);
	//
	// object that is derived from cppad_mixed
	mixed_derived mixed_object( n_fixed, n_random, y, sigma_u, sigma_y );
	mixed_object.initialize(fixed_vec, random_vec);
	//
	// lower and upper limits for random effects
	vector<double> random_lower(n_random), random_upper(n_random);
	for(size_t i = 0; i < n_random; i++)
	{	random_lower[i] = -inf;
		random_upper[i] = +inf;
	}
	//
	// optimize the random effects
	std::string options;
	options += "Integer print_level 0\n";
	options += "String  sb          yes\n";
	options += "String  derivative_test second-order\n";
	vector<double> uhat(n_random);
	uhat = mixed_object.optimize_random(
			options, fixed_vec, random_lower, random_upper, random_vec
	);
	//
	// compute the derivative of the random part of objective
	vector<double> r_fixed(n_fixed);
	mixed_object.ranobj_grad(fixed_vec, uhat, r_fixed);
	//
	// check the derivative of the random part of objective
	uhat = mixed_object.optimize_random(
			options, fixed_vec, random_lower, random_upper, random_vec
	);
	double r_theta = mixed_object.r_theta(fixed_vec[0], uhat[0]);
	//
	ok &= CppAD::abs( r_fixed[0] / r_theta - 1.0 ) < eps;
	//
	return ok;
}
// END C++
