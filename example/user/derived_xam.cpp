// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin derived_xam.cpp$$
$spell
	CppAD
	cppad
	interp
	xam
$$

$section mixed_cppad Derived Class: Example and Test$$

$code
$verbatim%example/user/derived_xam.cpp
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
		// implementation of ran_likelihood
		template <class Float>
		vector<Float> implement_ran_likelihood(
			const vector<Float>& theta  ,
			const vector<Float>& u      )
		{	vector<Float> vec(1);

			// compute this factor once
			Float sqrt_2pi = Float( CppAD::sqrt( 8.0 * CppAD::atan(1.0) ) );

			// initialize summation
			vec[0] = Float(0.0);

			// for each data and random effect
			for(size_t i = 0; i < y_.size(); i++)
			{	Float mu     = u[i];
				Float sigma  = theta[i];
				Float res    = (y_[i] - mu) / sigma;

				// This is a Gaussian term, so entire density is smooth
				vec[0]  += log(sqrt_2pi * sigma) + res * res / Float(2.0);
			}
			return vec;
		}
		// implementation of fix_likelihood
		template <class Float>
		vector<Float> implement_fix_likelihood(
			const vector<Float>& theta  )
		{	vector<Float> vec(1);

			// initialize part of log-density that is smooth
			vec[0] = Float(0.0);

			// compute these factors once
			Float one    = Float(1.0);
			Float sqrt_2 = CppAD::sqrt( Float(2.0) );

			for(size_t i = 0; i < y_.size(); i++)
			{	Float sigma  = theta[i];

				// This is a Laplace term, so using Laso for this variable
				vec[0] += CppAD::log( sqrt_2 * sigma);

				// part of the density that needs absolute value
				vec.push_back(sqrt_2 * (sigma - one) );
			}
			return vec;
		}
		// ------------------------------------------------------------------
		virtual vector<a2_double> ran_likelihood(
			const vector<a2_double>& fixed_vec  ,
			const vector<a2_double>& random_vec )
		{	return implement_ran_likelihood(fixed_vec, random_vec); }
		virtual vector<a1_double> ran_likelihood(
			const vector<a1_double>& fixed_vec  ,
			const vector<a1_double>& random_vec )
		{	return implement_ran_likelihood(fixed_vec, random_vec); }
		//
		virtual vector<a1_double> fix_likelihood(
			const vector<a1_double>& fixed_vec  )
		{	return implement_fix_likelihood(fixed_vec); }
		//
		virtual void fatal_error(const std::string& error_message)
		{	std::cerr << "fatal_error = " << error_message << std::endl;
			std::exit(1);
		}
		//
		virtual void warning(const std::string& warning_message)
		{	std::cerr << "warning = " << warning_message << std::endl;
		}
	};
}

bool derived_xam(void)
{
	bool   ok  = true;
	double pi  = 4.0 * std::atan(1.0);
	double eps = 100. * std::numeric_limits<double>::epsilon();

	typedef cppad_mixed::a1_double a1_double;
	typedef cppad_mixed::a2_double a2_double;

	size_t n_data   = 10;
	size_t n_fixed  = n_data;
	size_t n_random = n_data;
	vector<double>    data(n_data);
	vector<double>    fixed_vec(n_fixed), random_vec(n_random);
	vector<a1_double> a1_fixed(n_fixed), a1_random(n_random);
	vector<a2_double> a2_fixed(n_fixed), a2_random(n_random);

	for(size_t i = 0; i < n_data; i++)
	{	data[i]       = double(i + 1);
		//
		fixed_vec[i]  = 1.5;
		a1_fixed[i]   = a1_double( fixed_vec[i] );
		a2_fixed[i]   = a2_double( fixed_vec[i] );
		//
		random_vec[i] = 0.0;
		a1_random[i]  = a1_double( random_vec[i] );
		a2_random[i]  = a2_double( random_vec[i] );
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed = true;
	mixed_derived mixed_object(n_fixed, n_random, quasi_fixed, data);
	mixed_object.initialize(fixed_vec, random_vec);

	// Evaluate the random likelihood
	vector<a2_double> a2_vec(1);
	a2_vec = mixed_object.implement_ran_likelihood(a2_fixed, a2_random);

	// check the random likelihood
	double sum = 0.0;
	for(size_t i = 0; i < n_data; i++)
	{	double mu     = random_vec[i];
		double sigma  = fixed_vec[i];
		double res    = (data[i] - mu) / sigma;
		sum          += (std::log(2 * pi * sigma * sigma) + res * res) / 2.0;
	}
	ok &= abs( a2_vec[0] / a2_double(sum) - a2_double(1.0) ) < eps;

	// Evaluate the fixed likelihood
	vector<a1_double> a1_vec(1 + n_fixed);
	a1_vec = mixed_object.fix_likelihood(a1_fixed);

	// check the fixed likelihood
	sum = 0.0;
	for(size_t j = 0; j < n_fixed; j++)
	{	double sigma  = fixed_vec[j];
		sum          += std::log( std::sqrt(2.0) * sigma );
		//
		double check  = std::sqrt(2.0) * (sigma - 1.0);
		ok  &= abs( a1_vec[1+j] / check - 1.0 ) < eps;
	}
	ok &= abs( a1_vec[0] / sum - 1.0 ) < eps;

	return ok;
}
// END C++
