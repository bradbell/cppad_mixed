// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin optimize_random_xam.cpp$$
$spell
	CppAD
	cppad
	interp
	xam
$$

$section Optimize Random Effects: Example and Test$$

$code
$verbatim%example/user/optimize_random_xam.cpp
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

	class mixed_derived : public cppad_mixed {
	private:
		const vector<double>& y_;
	public:
		// constructor
		mixed_derived(
			size_t n_fixed                    ,
			size_t n_random                   ,
			bool   quasi_fixed                ,
			const vector<double>& y           )
			:
			cppad_mixed(n_fixed, n_random, quasi_fixed) ,
			y_(y)
		{ }
	private:
		// implementation of ran_like
		template <class Float>
		vector<Float> implement_ran_like(
			const vector<Float>& theta  ,
			const vector<Float>& u      )
		{	vector<Float> vec(1);

			// initialize part of log-density that is always smooth
			vec[0] = Float(0.0);

			// compute factor once
			Float sqrt_2pi = Float( CppAD::sqrt( 8.0 * CppAD::atan(1.0) ) );

			for(size_t i = 0; i < y_.size(); i++)
			{	Float mu     = u[i];
				Float sigma  = theta[i];
				Float res    = (y_[i] - mu) / sigma;

				// Gaussian likelihood
				vec[0]  += (sqrt_2pi * log(sigma) + res*res) / Float(2.0);
			}
			return vec;
		}
	public:
		// -------------------------------------------------------------------
		// User defined virtual functions
		virtual vector<a2_double> ran_like(
			const vector<a2_double>& fixed_vec  ,
			const vector<a2_double>& random_vec )
		{	return implement_ran_like(fixed_vec, random_vec); }
		virtual vector<a1_double> ran_like(
			const vector<a1_double>& fixed_vec  ,
			const vector<a1_double>& random_vec )
		{	return implement_ran_like(fixed_vec, random_vec); }
		//
		// improper constant prior
		virtual vector<a1_double> fix_like(
			const vector<a1_double>& fixed_vec  )
		{	a1d_vector vec(1);
			vec[0] = 0.0;
			return vec;
		}
		//
		virtual vector<a1_double> fix_constraint(
			const vector<a1_double>& fixed_vec  )
		{	return a1d_vector(0); } // empty vector
		//
		virtual void fatal_error(const std::string& error_message)
		{	std::cerr << "Error: " << error_message << std::endl;
			std::exit(1);
		}
		//
		virtual void warning(const std::string& warning_message)
		{	std::cerr << "Warning: " << warning_message << std::endl;
		}
		// ------------------------------------------------------------------
	};
}

bool optimize_random_xam(void)
{
	bool   ok = true;

	size_t n_data = 10;
	vector<double> data(n_data), fixed_vec(n_data), random_in(n_data);

	for(size_t i = 0; i < n_data; i++)
	{	data[i]      = double(i + 1);
		fixed_vec[i] = 1.0;
		random_in[i] = 0.0;
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed = true;
	mixed_derived mixed_object(n_data, n_data, quasi_fixed, data);
	mixed_object.initialize(fixed_vec, random_in);

	// lower and upper limits for random effects
	double inf = std::numeric_limits<double>::infinity();
	vector<double> random_lower(n_data), random_upper(n_data);
	for(size_t i = 0; i < n_data; i++)
	{	random_lower[i] = -inf;
		random_upper[i] = +inf;
	}

	// determine the optimal random effects
	std::string options;
	options += "Integer print_level 0\n";
	options += "String  sb          yes\n";
	options += "String  derivative_test second-order\n";
	vector<double> random_out = mixed_object.optimize_random(
		options, fixed_vec, random_lower, random_upper, random_in
	);

	// check the result
	for(size_t i = 0; i < n_data; i++)
	{	// debugging print out
		// std::cout << random_out[i] / data[i] - 1.0 << std::endl;
		ok &= CppAD::abs(random_out[i] / data[i] - 1.0) < 1e-10;
	}

	return ok;
}
// END C++
