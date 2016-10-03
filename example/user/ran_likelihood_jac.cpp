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
$begin ran_likelihood_jac.cpp$$
$spell
	CppAD
	cppad
	interp
	xam
	Jacobian
$$

$section Random Likelihood Jacobian: Example and Test$$

$code
$srcfile%example/user/ran_likelihood_jac.cpp
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
		// ------------------------------------------------------------------
		// implementation of ran_likelihood
		virtual vector<a2_double> ran_likelihood(
			const vector<a2_double>& theta  ,
			const vector<a2_double>& u      )
		{	vector<a2_double> vec(1);

			// sqrt_2pi = CppAD::sqrt( 8.0 * CppAD::atan(1.0) );

			// initialize summation
			vec[0] = a2_double(0.0);

			// for each data and random effect
			for(size_t i = 0; i < y_.size(); i++)
			{	a2_double mu     = u[i];
				a2_double sigma  = theta[i];
				a2_double res    = (y_[i] - mu) / sigma;

				// This is a Gaussian term, so entire density is smooth
				vec[0]  += log(sigma) + res * res / a2_double(2.0);
				// following term does not depend on fixed or random effects
				// vec[0]  += log(sqrt_2pi);
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
	vector<a2_double> a2_fixed(n_fixed), a2_random(n_random);

	for(size_t i = 0; i < n_data; i++)
	{	data[i]       = double(i + 1);
		//
		fixed_vec[i]  = 1.5;
		a1_fixed[i]   = fixed_vec[i];
		a2_fixed[i]   = a1_fixed[i];
		//
		random_vec[i] = 0.0;
		a1_random[i]  = random_vec[i];
		a2_random[i]  = a1_random[i];
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed = true;
	CppAD::mixed::sparse_mat_info A_info; // empty matrix
	mixed_derived mixed_object(n_fixed, n_random, quasi_fixed, A_info, data);
	mixed_object.initialize(fixed_vec, random_vec);

	// record Evaluation of a1_double version of random likelihood
	CppAD::Independent(a2_random);
	vector<a2_double> a2_vec(1);
	a2_vec = mixed_object.ran_likelihood(a2_fixed, a2_random);
	CppAD::ADFun<a1_double> a1_f(a2_random, a2_vec);

	// Evaluate ran_likelihood_jac
	vector<a1_double> a1_jac(n_random);
	a1_jac = mixed_object.ran_likelihood_jac(a1_fixed, a1_random);

	// Compute the Jacobian using AD
	vector<a1_double> check = a1_f.Jacobian(a1_random);

	// check the jacobian values
	// (mixed_object.initilaize also does this check)
	for(size_t i = 0; i < n_data; i++)
		ok &= CppAD::NearEqual( Value(a1_jac[i]), Value(check[i]), eps, eps);

	return ok;
}
// END C++
