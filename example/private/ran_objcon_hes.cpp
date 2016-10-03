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
$begin ran_objcon_hes.cpp$$
$spell
	objcon
	CppAD
	ran_obj
	cppad
	obj
	hes
	eval
	interp
	xam
$$

$section ran_objcon_hes: Example and Test$$

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
It follows that the Laplace approximation is exact and
$latex \[
	\B{p}( y_i | \theta ) \sim \B{N} \left( \theta_0 , 1 + \theta_1^2 \right)
\] $$

$code
$srcfile%example/private/ran_objcon_hes.cpp
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

bool ran_objcon_hes_xam(void)
{
	bool   ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();
	double sqrt_2pi = CppAD::sqrt(8.0 * CppAD::atan(1.0) );

	size_t n_data   = 10;
	size_t n_fixed  = 2;
	size_t n_random = n_data;
	vector<double> data(n_data), fixed_vec(n_fixed), random_vec(n_random);
	vector<double> beta(n_fixed), theta(n_fixed), uhat(n_random);

	fixed_vec[0] = 2.0;
	fixed_vec[1] = 0.5;
	for(size_t i = 0; i < n_data; i++)
	{	data[i]       = double(i + 1);
		random_vec[i] = i / double(n_data);
	}

	// object that is derived from cppad_mixed
	CppAD::mixed::sparse_mat_info A_info; // empty matrix
	mixed_derived mixed_object(n_fixed, n_random, A_info, data);
	mixed_object.initialize(fixed_vec, random_vec);

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

	// compute Hessian of random part of Laplace approximation
	vector<double> weight(1); // no random constraints
	weight[0] = 1.0;
	vector<size_t> row, col;
	vector<double> val;
	mixed_object.ran_objcon_hes(fixed_vec, uhat, weight, row, col, val);

	// check size of result vectors
	size_t K = row.size();
	ok &= col.size() == K;
	ok &= val.size() == K;

	// For this case the Laplace approximation is exactly equal the integral
	// p(y | theta ) = integral of p(y | theta , u) p(u | theta) du
	// record the function L(theta) = p(y | theta)
	typedef AD<double> a1_double;
	vector<a1_double> a1_fixed_vec(n_fixed), a1_L(1);
	a1_fixed_vec[0]  = fixed_vec[0];
	a1_fixed_vec[1]  = fixed_vec[1];
	CppAD::Independent(a1_fixed_vec);
	a1_double mu     = a1_fixed_vec[0];
	a1_double sigma  = a1_fixed_vec[1];
	a1_double delta  = CppAD::sqrt( sigma * sigma + 1.0 );
	a1_L[0] = 0.0;
	for(size_t i = 0; i < n_data; i++)
	{	a1_double res  = (data[i] - mu) / delta;
		a1_L[0]       += log(sqrt_2pi * delta) + res*res / 2.0;
	}
	CppAD::ADFun<double> f(a1_fixed_vec, a1_L);

	// compute Hessian of L
	vector<double> hes = f.Hessian(fixed_vec, 0);

	// check results
	ok &= hes.size() == 4;
	size_t check_count = 0;
	size_t non_zero    = 0;
	for(size_t i = 0; i < 2; i++)
	{	for(size_t j = 0; j <= i; j++)
		{	// only check lower triangle non-zero values
			if( hes[ i * 2 + j ] != 0.0 )
			{	non_zero++;
				for(size_t k = 0; k < K; k++)
				{	if( row[k] == i && col[k] == j )
					{	ok &= abs( val[k] / hes[i*n_fixed+j] - 1.0 ) <= eps;
						check_count++;
					}
				}
			}
		}
	}
	// should be only two non
	ok &= check_count == non_zero;
	ok &= K == non_zero;

	return ok;
}
// END C++
