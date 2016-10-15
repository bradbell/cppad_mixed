// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin ran_obj_jac.cpp$$
$spell
	jac
	CppAD
	ran_obj
	cppad
	obj
	interp
	xam
$$

$section ran_obj_jac: Example and Test$$

$head Private$$
This example is not part of the
$cref/cppad_mixed public API/public/$$.

$head Model$$
$latex \[
\B{p}( y_i | \theta , u ) \sim \B{N} ( \theta_0 + \theta_1 u_i, \theta_2^2 )
\] $$
$latex \[
	\B{p}( u_i | \theta ) \sim \B{N} ( 0 , 1 )
\] $$
It follows that the Laplace approximation is exact and
$latex \[
\B{p}( y_i | \theta )
\sim
\B{N} \left( \theta_0 , \theta_1^2 + \theta_2^2 \right)
\] $$

$code
$srcfile%example/private/ran_obj_jac.cpp
	%0%// BEGIN C++%// END C++%1%$$
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
	using CppAD::mixed::sparse_mat_info;
	//
	typedef AD<double>    a1_double;
	typedef AD<a1_double> a2_double;

	class mixed_derived : public cppad_mixed {
	private:
		const vector<double>& y_;
	public:
		// constructor
		mixed_derived(
			size_t n_fixed                    ,
			size_t n_random                   ,
			const vector<double>& y           )
			:
			// quasi_fixed = false
			cppad_mixed(n_fixed, n_random, false) ,
			y_(y)
		{	assert( n_fixed == 3);
		}
		// implementation of ran_likelihood
		virtual vector<a2_double> ran_likelihood(
			const vector<a2_double>& theta  ,
			const vector<a2_double>& u      )
		{	vector<a2_double> vec(1);

			// initialize part of log-density that is always smooth
			vec[0] = a2_double(0.0);

			// pi
			a2_double sqrt_2pi = a2_double(
				 CppAD::sqrt(8.0 * CppAD::atan(1.0)
			));

			for(size_t i = 0; i < y_.size(); i++)
			{	a2_double mu     = theta[0] + theta[1] * u[i];
				a2_double sigma  = theta[2];
				a2_double res    = (y_[i] - mu) / sigma;

				// p(y_i | u, theta)
				vec[0] += log(sqrt_2pi * sigma) + res*res / a2_double(2.0);

				// p(u_i | theta)
				vec[0] += log(sqrt_2pi) + u[i] * u[i] / a2_double(2.0);
			}
			return vec;
		}
	public:
		//
	};
}

bool ran_obj_jac_xam(void)
{
	bool   ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();

	size_t n_data   = 10;
	size_t n_fixed  = 3;
	size_t n_random = n_data;
	vector<double> data(n_data), fixed_vec(n_fixed), random_vec(n_random);
	vector<double> uhat(n_random);

	fixed_vec[0] = 2.0;
	fixed_vec[1] = 0.5;
	fixed_vec[2] = 1.0;
	for(size_t i = 0; i < n_data; i++)
	{	data[i]       = double(i + 1);
		random_vec[i] = i / double(n_data);
	}

	// object that is derived from cppad_mixed
	CppAD::mixed::sparse_mat_info A_info; // empty matrix
	mixed_derived mixed_object(n_fixed, n_random, data);
	mixed_object.initialize(fixed_vec, random_vec, A_info);

	// lower and upper limits for random effects
	double inf = std::numeric_limits<double>::infinity();
	vector<double> random_lower(n_random), random_upper(n_random);
	for(size_t i = 0; i < n_random; i++)
	{	random_lower[i] = -inf;
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
	vector<double> r_fixed(n_fixed);
	mixed_object.ran_obj_jac(fixed_vec, uhat, r_fixed);

	// For this case the Laplace approximation is exactly equal the integral
	// p(y | theta ) = integral of p(y | theta , u) p(u | theta) du
	// Furthermore p(y | theta ) is simple to calculate directly
	double mu         = fixed_vec[0];
	double alpha      = fixed_vec[1];
	double gamma      = fixed_vec[2];
	//
	double d_mu_0     = 1.0;
	double d_alpha_1  = 1.0;
	double d_gamma_2  = 1.0;
	//
	double delta      = CppAD::sqrt( alpha * alpha + gamma * gamma );
	double d_delta_1  = alpha * d_alpha_1 / delta;
	double d_delta_2  = gamma * d_gamma_2 / delta;
	//
	// double sum     = 0.0;
	double d_sum_0    = 0.0;
	double d_sum_1    = 0.0;
	double d_sum_2    = 0.0;
	for(size_t i = 0; i < n_data; i++)
	{	double res        = (data[i] - mu) / delta;
		double d_res_0    = - d_mu_0 / delta;
		double d_res_1    = - d_delta_1 * (data[i] - mu) / ( delta * delta) ;
		double d_res_2    = - d_delta_2 * (data[i] - mu) / ( delta * delta) ;
		//
		// double square  = res * res / 2.0;
		double d_square_0 = d_res_0 * res;
		double d_square_1 = d_res_1 * res;
		double d_square_2 = d_res_2 * res;
		//
		// double logdelta = CppAD::log(delta);
		double d_log_1     = d_delta_1 / delta;
		double d_log_2     = d_delta_2 / delta;
		//
		// sum           += logdelta + square;
		d_sum_0          += d_square_0;
		d_sum_1          += d_log_1 + d_square_1;
		d_sum_2          += d_log_2 + d_square_2;
	}
	ok &= abs( r_fixed[0] / d_sum_0 - 1.0 )  < eps;
	ok &= abs( r_fixed[1] / d_sum_1 - 1.0 )  < eps;
	ok &= abs( r_fixed[1] / d_sum_2 - 1.0 )  < eps;

	return ok;
}
// END C++
