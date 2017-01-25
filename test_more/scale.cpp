// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
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
		size_t                n_fixed_;
		const vector<double>& z_;
	public:
		// constructor
		mixed_derived(
			size_t n_fixed                    ,
			size_t n_random                   ,
			bool   quasi_fixed                ,
			const vector<double>& z           ) :
			cppad_mixed(n_fixed, n_random, quasi_fixed) ,
			n_fixed_(n_fixed)                                      ,
			z_(z)
		{	assert(z.size() == n_fixed); }
		// implementation of fix_likelihood as p(z|theta) * p(theta)
		virtual vector<a1_double> fix_likelihood(
			const vector<a1_double>& fixed_vec  )
		{
			// initialize log-density
			vector<a1_double> vec(1);
			vec[0] = a1_double(0.0);

			for(size_t j = 0; j < n_fixed_; j++)
			{	// case where partial w.r.t. theta_j does not exist
				// when theta_j == z_j
				a1_double sum_var  = 0.0;
				double    sum_data = 0.0;
				for(size_t k = 0; k <= j; k++)
				{	sum_var  += fixed_vec[k];
					sum_data += z_[k];
				}
				a1_double res = 1e+8 * (sum_var - sum_data);
				vec[0] += res * res / 2.0;
			}
			return vec;
		}
	};
}

bool scale(void)
{
	bool   ok = true;
	double inf = std::numeric_limits<double>::infinity();
	double tol = 1e-8;
	//
	// data
	size_t n_fixed  = 5;
	vector<double> z(n_fixed);
	double scale = 1.0;
	for(size_t i = 0; i < n_fixed; i++)
		z[i] = scale * double(i+1);

	// fixed effects
	vector<double>
		fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
	for(size_t j = 0; j < n_fixed; j++)
	{	fixed_lower[j] = - inf;
		fixed_in[j]    = 0.0;
		fixed_upper[j] = inf;
	}
	//
	// no random effects
	size_t n_random = 0;
	vector<double> random_in(0);
	//
	// no constriants
	vector<double> fix_constraint_lower(0), fix_constraint_upper(0);

	// object that is derived from cppad_mixed
	bool quasi_fixed = true;
	CppAD::mixed::sparse_mat_info A_info; // empty matrix
	mixed_derived mixed_object(n_fixed, n_random, quasi_fixed, z);
	mixed_object.initialize(fixed_in, random_in, A_info);

	// optimize the fixed effects using quasi-Newton method
	std::string fixed_ipopt_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           first-order\n"
		"String  derivative_test_print_all yes\n"
		"Numeric tol                       1e-8\n"
	;
	std::string random_ipopt_options =
		"Integer print_level 0\n"
		"String  sb          yes\n"
		"String  derivative_test second-order\n"
	;
	vector<double> random_lower(n_random), random_upper(n_random);
	for(size_t i = 0; i < n_random; i++)
	{	random_lower[i] = -inf;
		random_upper[i] = +inf;
	}
	CppAD::mixed::fixed_solution solution = mixed_object.optimize_fixed(
		fixed_ipopt_options,
		random_ipopt_options,
		fixed_lower,
		fixed_upper,
		fix_constraint_lower,
		fix_constraint_upper,
		fixed_in,
		random_lower,
		random_upper,
		random_in
	);
	vector<double> fixed_out = solution.fixed_opt;
	//
	for(size_t j = 0; j < n_fixed; j++)
	{	double err = fabs( fixed_out[j] / z[j] - 1.0 );
		ok &= err <= tol;
	}

	return ok;
}
// END C++
