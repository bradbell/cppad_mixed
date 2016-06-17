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
$begin optimize_random_xam.cpp$$
$spell
	CppAD
	cppad
	interp
	xam
$$

$section Optimize Random Effects: Example and Test$$

$code
$srcfile%example/user/optimize_random_xam.cpp
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
			bool   quasi_fixed                ,
			const  sparse_mat_info& A_info    ,
			const vector<double>& y           )
			:
			cppad_mixed(n_fixed, n_random, quasi_fixed, A_info) ,
			y_(y)
		{ }
		// implementation of ran_likelihood
		virtual vector<a2_double> ran_likelihood(
			const vector<a2_double>& theta  ,
			const vector<a2_double>& u      )
		{	vector<a2_double> vec(1);

			// initialize part of log-density that is always smooth
			vec[0] = a2_double(0.0);

			// sqrt_2pi = CppAD::sqrt( 8.0 * CppAD::atan(1.0) );

			for(size_t i = 0; i < y_.size(); i++)
			{	a2_double mu     = u[i];
				a2_double sigma  = theta[i];
				a2_double res    = (y_[i] - mu) / sigma;

				// Gaussian likelihood
				vec[0]  += log(sigma) + res * res / a2_double(2.0);
				// following term does not depend on fixed or random effects
				// vec[0]  += log(sqrt_2pi);
			}
			return vec;
		}
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
	CppAD::mixed::sparse_mat_info A_info; // empty matrix
	mixed_derived mixed_object(n_data, n_data, quasi_fixed, A_info, data);
	mixed_object.initialize(fixed_vec, random_in);

	// lower and upper limits for random effects
	double inf = std::numeric_limits<double>::infinity();
	vector<double> random_lower(n_data), random_upper(n_data);
	for(size_t i = 0; i < n_data; i++)
	{	random_lower[i] = -inf;
		random_upper[i] = +inf;
	}

	// -----------------------------------------------------------------------
	// use ipopt to determine the optimal random effects
	std::string ipopt_options;
	ipopt_options += "Integer print_level 0\n";
	ipopt_options += "String  sb          yes\n";
	ipopt_options += "String  derivative_test second-order\n";
	vector<double> random_out = mixed_object.optimize_random(
		ipopt_options, fixed_vec, random_lower, random_upper, random_in
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
