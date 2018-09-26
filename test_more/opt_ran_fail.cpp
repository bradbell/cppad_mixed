// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>

// optimizing random effects fails when fixed effects close to solution
namespace {
	using CppAD::vector;
	using CppAD::log;
	using CppAD::AD;
	//
	using CppAD::mixed::d_sparse_rcv;
	using CppAD::mixed::d_vector;
	//
	double small = 1e-5;
	//
	class mixed_derived : public cppad_mixed {
	private:
		const d_vector&       y_;
	public:
		// constructor
		mixed_derived(
			size_t                 n_fixed       ,
			size_t                 n_random      ,
			bool                   quasi_fixed   ,
			bool                   bool_sparsity ,
			const d_sparse_rcv&    A_rcv         ,
			const d_vector&        y             ) :
			cppad_mixed(
				n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
			),
			y_(y)
		{ }
		// implementation of ran_likelihood
		template <typename Vector>
		Vector template_ran_likelihood(
			const Vector&         theta  ,
			const Vector&         u      )
		{	typedef typename Vector::value_type scalar;

			Vector vec(1);

			// sum of residual squared
			vec[0]            = scalar(0.0);
			scalar sum_sq  = scalar(0.0);
			for(size_t i = 0; i < y_.size(); i++)
			{	// data term
				scalar mu     = theta[i] + u[i];
				scalar res    = y_[i] - mu;

				// p(y|theta,u)
				vec[0] += res * res;

				// add to sum (how far fixed effects from solution)
				res     = y_[i] - theta[i];
				sum_sq += res * res;

				// p(u|theta)
				vec[0] += u[i] * u[i];
			}

			// return nan when sum of squares is less than 1e-4
			scalar vec_0    = CppAD::numeric_limits<scalar>::quiet_NaN();
			scalar a2_small = scalar( small );
			vec_0  = CppAD::CondExpGt(sum_sq, + a2_small, vec[0], vec_0);
			vec[0] = vec_0;
			//
			return vec;
		}
		// a1_vector version of ran_likelihood
		virtual a1_vector ran_likelihood(
			const a1_vector& fixed_vec, const a1_vector& random_vec
		)
		{	return template_ran_likelihood( fixed_vec, random_vec ); }
		// we expect to get a warnings
		virtual void warning(const std::string& warning_message)
		{ }
	};
}

bool opt_ran_fail(void)
{
	bool   ok  = true;
	double inf = std::numeric_limits<double>::infinity();

	size_t n_data = 5;
	d_vector data(n_data);
	d_vector random_lower(n_data), random_in(n_data), random_upper(n_data);
	d_vector fixed_lower(n_data), fixed_in(n_data), fixed_upper(n_data);
	for(size_t i = 0; i < n_data; i++)
	{	data[i]      = double(i + 1);
		//
		fixed_lower[i]  = -inf;
		fixed_in[i]     = 0.0;
		fixed_upper[i]  = +inf;
		//
		random_lower[i] = -inf;
		random_in[i]    = data[i];
		random_upper[i] = +inf;
	}
	d_vector fix_constraint_lower(0), fix_constraint_upper(0);

	// object that is derived from cppad_mixed
	bool quasi_fixed   = true;
	bool bool_sparsity = true;
	d_sparse_rcv A_rcv; // empty matrix
	size_t n_fixed  = n_data;
	size_t n_random = n_data;
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, data
	);
	mixed_object.initialize(fixed_in, random_in);

	// optimize the fixed effects using quasi-Newton method
	std::string fixed_ipopt_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           adaptive\n"
		"String  derivative_test_print_all yes\n"
		"Integer max_iter                  15\n"
		"Numeric tol                       1e-8\n"
	;
	std::string random_ipopt_options =
		"Integer print_level     0\n"
		"String  sb              yes\n"
		"String  derivative_test second-order\n"
		"Numeric tol             1e-8\n"
	;
	d_vector fixed_scale = fixed_in;
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
	d_vector fixed_out = solution.fixed_opt;

	// check that residual got to small before aborting
	double sum_sq = 0.0;
	for(size_t i = 0; i < n_data; i++)
	{	double res = data[i] - fixed_out[i];
		sum_sq    += res * res;
	}
	ok &= sum_sq <= 2.0 * small;
	//
	return ok;
}
// END C++
