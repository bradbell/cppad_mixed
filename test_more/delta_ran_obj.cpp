// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
// Finite difference test of ran_obj_jac.cpp
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>

namespace {
	using CppAD::vector;
	using CppAD::log;
	using CppAD::exp;
	using CppAD::AD;
	using CppAD::mixed::d_sparse_rcv;
	//
	using CppAD::mixed::a1_double;

	class mixed_derived : public cppad_mixed {
	private:
		const double y_;
		const double sigma_u_, sigma_y_;
	public:
		// constructor
		mixed_derived(
			size_t n_fixed         ,
			size_t n_random        ,
			const d_sparse_rcv A_rcv  ,
			double y               ,
			double sigma_u         ,
			double sigma_y         ) :
			// quasi_fixed=false, bool_sparsity=true,
			cppad_mixed(n_fixed, n_random, false, true, A_rcv) ,
			y_(y) , sigma_u_(sigma_u), sigma_y_(sigma_y)
		{	assert( n_fixed == 1 );
			assert( n_random == 2 );
		}
	public:
		// ------------------------------------------------------------------
		// implementation of ran_likelihood as p(y|theta, u) * p(u|theta)
		template <typename Vector>
		Vector template_ran_likelihood(
			const Vector& fixed_vec  ,
			const Vector& random_vec )
		{	typedef typename Vector::value_type scalar;

			scalar theta = fixed_vec[0];
			scalar sum_u = random_vec[0] + random_vec[1];

			// initialize log-density
			Vector vec(1);
			vec[0] = scalar(0.0);

			// compute this factors once
			scalar sqrt_2pi = scalar(
				 CppAD::sqrt( 8.0 * CppAD::atan(1.0)
			));

			// Data term
			scalar res  = (y_ - exp(sum_u) * theta) / sigma_y_;
			vec[0]    += log(sqrt_2pi * sigma_y_ ) + res * res / scalar(2.0);

			// prior for u
			for(size_t i = 0; i < 2; i++)
			{	res     = random_vec[i] / sigma_u_;
				vec[0] += log(sqrt_2pi * sigma_u_ ) + res * res / scalar(2.0);
			}
			return vec;
		}
		// a1_vector version of ran_likelihood
		virtual a1_vector ran_likelihood(
			const a1_vector& fixed_vec, const a1_vector& random_vec
		)
		{	return template_ran_likelihood( fixed_vec, random_vec ); }
		//
	};

}

bool delta_ran_obj(void)
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
	vector<double> fix_constraint_lower(0), fix_constraint_upper(0);
	//
	// object that is derived from cppad_mixed
	CppAD::mixed::d_sparse_rcv A_rcv; // empty matrix
	mixed_derived mixed_object(
			n_fixed, n_random, A_rcv, y, sigma_u, sigma_y
	);
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
	// factor f_{u,u} (theta, uhat)
	mixed_object.update_factor(fixed_vec, uhat);
	//
	// compute the derivative of the random part of objective
	vector<double> r_fixed(n_fixed);
	mixed_object.ran_obj_jac(fixed_vec, uhat, r_fixed);
	//
	// check using finite central differences
	for(size_t j = 0; j < n_fixed; j++)
	{	double theta_j = fixed_vec[j];
		fixed_vec[j]   = theta_j + 2.0 * eps;
		uhat = mixed_object.optimize_random(
			options, fixed_vec, random_lower, random_upper, random_vec
		);
		//
		// factor f_{u,u} (theta, uhat)
		mixed_object.update_factor(fixed_vec, uhat);
		//
		double r_plus  = mixed_object.ran_obj_eval(fixed_vec, uhat);
		fixed_vec[j]   = theta_j - 2.0 * eps;
		uhat = mixed_object.optimize_random(
			options, fixed_vec, random_lower, random_upper, random_vec
		);
		//
		// factor f_{u,u} (theta, uhat)
		mixed_object.update_factor(fixed_vec, uhat);
		//
		double r_minus = mixed_object.ran_obj_eval(fixed_vec, uhat);
		fixed_vec[j]   = theta_j;
		//
		double check = (r_plus - r_minus) / (4.0 * eps);
		//
		ok &= fabs( check / r_fixed[j] - 1.0) <= eps;
	}
	return ok;
}
// END C++
