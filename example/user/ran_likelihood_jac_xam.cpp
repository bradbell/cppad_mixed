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
$begin ran_likelihood_jac_xam.cpp$$
$spell
	CppAD
	cppad
	interp
	xam
$$

$section mixed_cppad Derived Class: Example and Test$$

$code
$srcfile%example/user/ran_likelihood_jac_xam.cpp
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
		// ------------------------------------------------------------------
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
		// ------------------------------------------------------------------
		// ran_likelihood_jac
		vector<a1_double> ran_likelihood_jac(
			const vector<a1_double>& theta  ,
			const vector<a1_double>& u      )
		{
			// return value
			vector<a1_double> vec(y_.size());

			// for each data and random effect
			for(size_t i = 0; i < y_.size(); i++)
			{	a1_double mu     = u[i];
				a1_double sigma  = theta[i];
				a1_double res    = (y_[i] - mu) / sigma;
				a1_double res_ui = - 1.0 / sigma;

				// This is a Gaussian term, so entire density is smooth
				vec[i]  =  res * res_ui;
			}
			return vec;
		}
		// ------------------------------------------------------------------
		// example a2 version of ran_likelihood
		virtual vector<a2_double> ran_likelihood(
			const vector<a2_double>& fixed_vec  ,
			const vector<a2_double>& random_vec )
		{	return implement_ran_likelihood(fixed_vec, random_vec); }
		// ------------------------------------------------------------------
		// example a1 version of ran_likelihood
		virtual vector<a1_double> ran_likelihood(
			const vector<a1_double>& fixed_vec  ,
			const vector<a1_double>& random_vec )
		{	return implement_ran_likelihood(fixed_vec, random_vec); }
	};
}

bool ran_likelihood_jac_xam(void)
{
	bool   ok  = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();
	typedef cppad_mixed::a1_double a1_double;

	size_t n_data   = 10;
	size_t n_fixed  = n_data;
	size_t n_random = n_data;
	vector<double>    data(n_data);
	vector<double>    fixed_vec(n_fixed), random_vec(n_random);
	vector<a1_double> a1_fixed(n_fixed), a1_random(n_random);

	for(size_t i = 0; i < n_data; i++)
	{	data[i]       = double(i + 1);
		//
		fixed_vec[i]  = 1.5;
		a1_fixed[i]   = a1_double( fixed_vec[i] );
		//
		random_vec[i] = 0.0;
		a1_random[i]  = a1_double( random_vec[i] );
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed = true;
	CppAD::mixed::sparse_mat_info A_info; // empty matrix
	mixed_derived mixed_object(n_fixed, n_random, quasi_fixed, A_info, data);
	mixed_object.initialize(fixed_vec, random_vec);

	// record Evaluation of a1_double version of random likelihood
	CppAD::Independent(a1_random);
	vector<a1_double> a1_fun(1);
	a1_fun = mixed_object.ran_likelihood(a1_fixed, a1_random);
	CppAD::ADFun<double> f(a1_random, a1_fun);

	// Evaluate ran_likelihood_jac
	vector<a1_double> a1_jac(n_random);
	a1_jac = mixed_object.ran_likelihood_jac(a1_fixed, a1_random);

	// Compute the Jacobian using AD
	vector<double> check = f.Jacobian(random_vec);

	// check the jacobian values
	// (when NDEBUG is not defined, cppad_mixed also does this check)
	for(size_t i = 0; i < n_data; i++)
		ok &= CppAD::NearEqual( Value(a1_jac[i]), check[i], eps, eps);

	return ok;
}
// END C++
