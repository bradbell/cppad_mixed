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
$begin logdet_jac.cpp$$
$spell
	jac
	CppAD
	cppad
	hes
	interp
	xam
	logdet
$$

$section logdet_jac: Example and Test$$


$head Private$$
This example is not part of the
$cref/cppad_mixed public API/public/$$.

$code
$srcfile%example/private/logdet_jac.cpp
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
			const  sparse_mat_info& A_info    ,
			const vector<double>& y           )
			:
			// quasi_fixed = false
			cppad_mixed(n_fixed, n_random, false, A_info) ,
			y_(y)
		{ }
		// implementation of ran_likelihood
		virtual vector<a2_double> ran_likelihood(
			const vector<a2_double>& theta  ,
			const vector<a2_double>& u      )
		{	vector<a2_double> vec(1);

			// initialize part of log-density that is always smooth
			vec[0] = a2_double(0.0);

			for(size_t i = 0; i < y_.size(); i++)
			{	a2_double mu     = u[i];
				a2_double sigma  = theta[i];
				a2_double res    = (y_[i] - mu) / sigma;

				// (do not need 2*pi inside of log)
				vec[0]  += (log(sigma) + res*res) / a2_double(2.0);
			}
			return vec;
		}
	public:
		//
		//
		virtual vector<a1_double> fix_likelihood(
			const vector<a1_double>& fixed_vec  )
		{	a1d_vector vec(1);
			vec[0] = 0.0;
			return vec;
		}
	};
}

bool logdet_jac_xam(void)
{
	bool   ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();

	size_t n_data   = 10;
	size_t n_fixed  = n_data;
	size_t n_random = n_data;
	vector<double> data(n_data);
	vector<double> theta(n_fixed), u(n_random);
	vector<double> fixed_vec(n_fixed), random_vec(n_random);

	for(size_t i = 0; i < n_data; i++)
	{	data[i]      = double( (i + 1) * (i + 1) );
		fixed_vec[i] = theta[i] = std::sqrt( double(i + 1) );
		random_vec[i] = u[i] = 0.0;
	}

	// object that is derived from cppad_mixed
	CppAD::mixed::sparse_mat_info A_info; // empty matrix
	mixed_derived mixed_object(n_fixed, n_random, A_info, data);
	mixed_object.initialize(theta, u);

	// factor f_{u,u} (thete, u)
	mixed_object.update_factor(fixed_vec, random_vec);

	// compute derivative of logdet of Hessian
	vector<double> logdet_fix(n_fixed), logdet_ran(n_random);
	mixed_object.logdet_jac(fixed_vec, random_vec, logdet_fix, logdet_ran);

	// Hessian_{i,j} = 1.0 / (theta[i] * theta[i]) if i == j
	//               = 0.0 otherwise
	// log( det( Hessian ) ) = - 2.0 * sum_i log( theta[i] )
	for(size_t i = 0; i < n_data; i++)
		ok    &= logdet_ran[i] == 0.0;
	for(size_t i = 0; i < n_data; i++)
	{	double check   = - 2.0  / fixed_vec[i];
		ok            &= abs( logdet_fix[i] / check - 1.0) <= eps;
	}

	return ok;
}
// END C++
