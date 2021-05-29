// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
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
	using CppAD::mixed::d_sparse_rcv;
	//
	using CppAD::mixed::a1_double;

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
			bool   bool_sparsity              ,
			const CppAD::mixed::d_sparse_rcv&    A_rcv,
			const vector<double>& z           ) :
			cppad_mixed(n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv),
			n_fixed_(n_fixed)                                      ,
			z_(z)
		{	assert(z.size() == n_fixed); }
		// implementation of fix_likelihood as p(z|theta) * p(theta)
		template <typename Vector>
		Vector template_fix_likelihood(
			const Vector& fixed_vec  )
		{	typedef typename Vector::value_type scalar;


			// initialize log-density
			Vector vec(1);
			vec[0] = scalar(0.0);

			for(size_t j = 0; j < n_fixed_; j++)
			{	scalar res = z_[j] - fixed_vec[j];
				vec[0] += exp( res * res );
			}
			//
			return vec;
		}
		// a1_vector version of fix_likelihood
		virtual a1_vector fix_likelihood(const a1_vector& fixed_vec)
		{	return template_fix_likelihood( fixed_vec ); }
	};
}

bool scale_two(void)
{
	bool   ok = true;
	double inf = std::numeric_limits<double>::infinity();
	double tol = 1e-10;
	//
	// data
	size_t n_fixed  = 5;
	vector<double> z(n_fixed);
	for(size_t i = 0; i < n_fixed; i++)
		z[i] = double(i+2);

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
	bool quasi_fixed  = false;
	bool bool_sparsity = false;
	CppAD::mixed::d_sparse_rcv A_rcv; // empty matrix
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, z
	);
	mixed_object.initialize(fixed_in, random_in);

	// optimize the fixed effects using quasi-Newton method
	std::string fixed_ipopt_options =
		"Integer print_level                5\n"
		"String  sb                         yes\n"
		"String  derivative_test            trace-adaptive\n"
		"String  derivative_test_print_all yes\n"
		"Numeric tol                        1e-10\n"
		"Integer max_iter                   50\n"
		"Integer limited_memory_max_history 10\n"
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
	vector<double> fixed_scale = fixed_in;
	CppAD::mixed::fixed_solution solution = mixed_object.optimize_fixed(
		fixed_ipopt_options,
		random_ipopt_options,
		fixed_lower,
		fixed_upper,
		fix_constraint_lower,
		fix_constraint_upper,
		fixed_scale,
		fixed_in,
		random_lower,
		random_upper,
		random_in
	);
	vector<double> fixed_out = solution.fixed_opt;
	//
	for(size_t j = 0; j < n_fixed; j++)
	{	double err = fabs( fixed_out[j] / z[j] - 1.0 );
		std::cout << "err[" << j << "]= " << err << "\n";
		ok &= err <= 5. * tol;
	}
	//
	return ok;
}
// END C++