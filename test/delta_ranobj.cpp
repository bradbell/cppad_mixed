// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
// Finite difference test of ranobj_grad.cpp
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad_mixed/cppad_mixed.hpp>

namespace {
	using CppAD::vector;
	using CppAD::log;
	using CppAD::exp;
	using CppAD::abs;
	using CppAD::AD;

	class mixed_derived : public cppad_mixed::cppad_mixed {
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
			// quasi_fixed = flase
			cppad_mixed::cppad_mixed(n_fixed, n_random, false) ,
			y_(y) , sigma_u_(sigma_u), sigma_y_(sigma_y)
		{	assert( n_fixed == 1 );
			assert( n_random == 2 );
		}
	public:
		// ------------------------------------------------------------------
		// implementation of ran_like as p(y|theta, u) * p(u|theta)
		template <class Float>
		vector<Float> implement_ran_like(
			const vector<Float>& fixed_vec  ,
			const vector<Float>& random_vec )
		{	Float theta = fixed_vec[0];
			Float sum_u = random_vec[0] + random_vec[1];

			// initialize log-density
			vector<Float> vec(1);
			vec[0] = Float(0.0);

			// compute this factors once
			Float sqrt_2pi = Float( CppAD::sqrt( 8.0 * CppAD::atan(1.0) ) );

			// Data term
			Float res  = (y_ - exp(sum_u) * theta) / sigma_y_;
			vec[0]    += log(sqrt_2pi * sigma_y_ ) + res * res / Float(2.0);

			// prior for u
			for(size_t i = 0; i < 2; i++)
			{	res     = random_vec[i] / sigma_u_;
				vec[0] += log(sqrt_2pi * sigma_u_ ) + res * res / Float(2.0);
			}
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
	};

}

bool delta_ranobj(void)
{
	bool   ok = true;
	double inf = std::numeric_limits<double>::infinity();
	double eps = sqrt( std::numeric_limits<double>::epsilon() );

	size_t n_fixed  = 1;
	size_t n_random = 2;
	double y        = 0.05;
	double sigma_u  = 0.2;
	double sigma_y  = 0.2 * y;
	vector<double> fixed_vec(n_fixed), random_vec(n_random);
	vector<double> fixed_lower(n_fixed), fixed_upper(n_fixed);
	fixed_lower[0] = -inf;
	fixed_vec[0]   = y;
	fixed_upper[0] = + inf;
	for(size_t i = 0; i < n_random; i++)
		random_vec[i] = 0.0;
	//
	// no constriants
	vector<double> constraint_lower(0), constraint_upper(0);
	//
	// object that is derived from cppad_mixed
	mixed_derived mixed_object( n_fixed, n_random, y, sigma_u, sigma_y );
	mixed_object.initialize(fixed_vec, random_vec);

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
	// check using finite central differences
	for(size_t j = 0; j < n_fixed; j++)
	{	double theta_j = fixed_vec[j];
		fixed_vec[j]   = theta_j + 2.0 * eps;
		uhat = mixed_object.optimize_random(
			options, fixed_vec, random_lower, random_upper, random_vec
		);
		double r_plus  = mixed_object.ranobj_eval(fixed_vec, uhat);
		fixed_vec[j]   = theta_j - 2.0 * eps;
		uhat = mixed_object.optimize_random(
			options, fixed_vec, random_lower, random_upper, random_vec
		);
		double r_minus = mixed_object.ranobj_eval(fixed_vec, uhat);
		fixed_vec[j]   = theta_j;
		//
		double check = (r_plus - r_minus) / (4.0 * eps);
		//
		ok &= CppAD::abs( check / r_fixed[j] - 1.0) <= eps;
	}
	return ok;
}
// END C++
