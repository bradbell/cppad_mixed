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
$begin ran_con_jac.cpp$$
$spell
	CppAD
	ran_con
	cppad
	xam
	jac
$$

$section ran_con_jac: Example and Test$$

$head Private$$
This example is not part of the
$cref/cppad_mixed public API/public/$$.

$head Model$$
$latex \[
	\B{p}( y_i | \theta , u ) \sim \B{N} ( u_i + \theta_0 , \theta_1^2 )
\] $$
$latex \[
	\B{p}( u_i | \theta ) \sim \B{N} ( 0 , 1 )
\] $$
$latex \[
	A = [ 1 , \cdots , 1 ]
\] $$
It follows that the $th i$$ component of the
optimal random effects minimize the function
$latex \[
	0.5 ( u_i + theta_0 - y_i )^2 / \theta_1^2 + 0.5 u_i^2
\] $$
Taking the derivative w.r.t $latex u_i$$ and setting it equal to zero
we have
$latex \[
	\hat{u}_i ( \theta ) = ( y_i - \theta_0 ) / ( 1 + \theta_1^2 )
\] $$
The random constraint function is
$latex \[
	\sum_i \hat{u}_i ( \theta )
	=
	\sum_i
	( y_i - \theta_0 ) / ( 1 + \theta_1^2 )
\] $$
The partial w.r.t. $latex \theta_0$$ is
$latex \[
	- \sum_i  1.0 / ( 1 + \theta_1^2 )
\] $$
The partial w.r.t. $latex \theta_1$$ is
$latex \[
	- 2 \sum_i ( y_i - \theta_0 ) \theta_1 / ( 1 + \theta_1^2 )
\] $$

$code
$srcfile%example/private/ran_con_jac.cpp
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
		{	assert( n_fixed == 2);
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
			{	a2_double mu     = u[i] + theta[0];
				a2_double sigma  = theta[1];
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

bool ran_con_jac_xam(void)
{
	bool   ok = true;

	size_t n_data   = 10;
	size_t n_fixed  = 2;
	size_t n_random = n_data;
	vector<double> data(n_data), fixed_vec(n_fixed), random_vec(n_random);
	vector<double> uhat(n_random);

	fixed_vec[0] = 2.0;
	fixed_vec[1] = 1.0;
	for(size_t i = 0; i < n_data; i++)
	{	data[i]       = double(i + 2);
		random_vec[i] = 0.0;
	}

	// lower and upper limits for random effects
	double inf = std::numeric_limits<double>::infinity();
	vector<double> random_lower(n_random), random_upper(n_random);
	for(size_t i = 0; i < n_random; i++)
	{	random_lower[i] = -inf;
		random_upper[i] = +inf;
	}

	// random constraint matrix
	CppAD::mixed::sparse_mat_info A_info;
	A_info.resize(n_random);
	for(size_t j = 0; j < n_random; j++)
	{	A_info.row[j] = 0;
		A_info.col[j] = j;
		A_info.val[j] = 1.0;
	}

	// object that is derived from cppad_mixed
	mixed_derived mixed_object(n_fixed, n_random, A_info, data);
	mixed_object.initialize(fixed_vec, random_vec);

	// optimize the random effects
	std::string options;
	options += "Integer print_level 0\n";
	options += "Numeric tol         1e-10\n";
	options += "String  sb          yes\n";
	options += "String  derivative_test second-order\n";
	uhat = mixed_object.optimize_random(
		options, fixed_vec, random_lower, random_upper, random_vec
	);

	// must factor f_{u,u} (theta, uhat)
	mixed_object.update_factor(fixed_vec, uhat);

	// compute sparstiy pattern for jacobian of random constraints
	CppAD::mixed::sparse_mat_info jac_info;
	mixed_object.ran_con_jac(fixed_vec, uhat, jac_info);

	// check number of possibly non_zero elements.
	ok &= jac_info.row.size() == 2;
	//
	// partial w.r.t. theta_0
	ok &= jac_info.row[0] == 0;
	ok &= jac_info.col[0] == 0;
	//
	// partial w.r.t. theta_1
	ok &= jac_info.row[1] == 0;
	ok &= jac_info.col[1] == 1;

	// Now compute the sparse Jacobian values
	mixed_object.ran_con_jac(fixed_vec, uhat, jac_info);
	//
	double jac_0 = jac_info.val[0];
	double jac_1 = jac_info.val[1];

	// check results
	double check_0 = 0.0;
	double check_1 = 0.0;
	for(size_t i = 0; i < n_random; i++)
	{	double theta_0   = fixed_vec[0];
		double theta_1   = fixed_vec[1];
		double theta_1sq = theta_1 * theta_1;
		check_0 -= 1.0 / (1.0 + theta_1sq );
		check_1 -= (data[i] - theta_0) * theta_1 / (1.0 + theta_1sq );
	}

	// used 1e-10 for random_optmize tolerance
	ok &= abs( jac_0 / check_0 - 1.0 ) < 1e-9;
	ok &= abs( jac_1 / check_1 - 1.0 ) < 1e-9;

	return ok;
}
// END C++
