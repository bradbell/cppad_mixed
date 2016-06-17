// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
// case with no cross terms; i.e, f_{u,theta} ( theta , u ) is zero

# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>

namespace {
	using CppAD::vector;
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
			{	a2_double mu     = u[i];
				a2_double sigma  = a2_double(0.5);
				a2_double res    = (y_[i] - mu) / sigma;

				// p(y_i | u, theta)
				vec[0] += CppAD::log(sqrt_2pi * sigma) + res*res / a2_double(2.0);

				// p(u_i | theta)
				vec[0] += CppAD::log(sqrt_2pi) + u[i] * u[i] / a2_double(2.0);
			}
			return vec;
		}
	public:
		//
	};
}

bool ran_obj_tst(void)
{
	bool   ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();

	size_t n_data   = 1;
	size_t n_fixed  = 2;
	size_t n_random = n_data;
	vector<double> data(n_data), fixed_vec(n_fixed), random_vec(n_random);
	vector<double> uhat(n_random);

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

	// factor f_{u,u} ( theta , uhat )
	mixed_object.update_factor(fixed_vec, uhat);

	// compute total derivative of random part of objective
	// vector<double> r_fixed(n_fixed);
	// mixed_object.ran_obj_jac(fixed_vec, uhat, r_fixed);


	// For this case the Laplace approximation is exactly equal the integral
	// p(y_i | theta ) = integral of p(y_i | theta , u) p(u | theta) du
	// Furthermore p(y_i | theta ) is N( theta[0], 1 + 0.5^2 )

	// check the random part of the objective
	double r        = mixed_object.ran_obj_eval(fixed_vec, uhat);
	double check    = 0.0;
	double sqrt_2pi = CppAD::sqrt(8.0 * CppAD::atan(1.0) );
	double sigma    = CppAD::sqrt( 1.0 + 0.5 * 0.5 );
	for(size_t i = 0; i < n_data; i++)
	{	double res   = data[i] / sigma;
		check       += CppAD::log(sqrt_2pi * sigma);
		check       += res * res / 2.0;
	}
	ok &= CppAD::abs( r / check - 1.0 ) <= eps;

	// check jacobian of random part of objective
	vector<double> r_fixed(n_fixed);
	mixed_object.ran_obj_jac(fixed_vec, uhat, r_fixed);
	//
	ok &= CppAD::abs( r_fixed[0] - 0.0 ) <= eps;
	ok &= CppAD::abs( r_fixed[1] - 0.0 ) <= eps;

	return ok;
}
