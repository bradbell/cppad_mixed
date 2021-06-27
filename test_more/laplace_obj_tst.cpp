/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
Test Laplace part of the the objective where non-linear and initial
random effects far from solution. (Used to faile before initialized
second order Lapalce objective at optimal random effects
(for initial fixed effects).
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/order2random.hpp>
# include <cppad/mixed/ran_like_hes.hpp>


namespace {
	using CppAD::vector;
	using CppAD::log;
	using CppAD::AD;
	using CppAD::mixed::d_sparse_rcv;
	using CppAD::mixed::d_vector;
	using CppAD::mixed::a1_double;
	using CppAD::mixed::a1_vector;
	using CppAD::mixed::sparse_rc;
	// -----------------------------------------------------------------------
	class mixed_derived : public cppad_mixed {
	private:
		const d_vector& y_;
	public:
		// constructor
		mixed_derived(
			size_t                               n_fixed       ,
			size_t                               n_random      ,
			bool                                 quasi_fixed   ,
			bool                                 bool_sparsity ,
			const CppAD::mixed::d_sparse_rcv&    A_rcv         ,
			const d_vector&                      y             ) :
			cppad_mixed(
				n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
			),
			y_(y)
		{	assert( n_fixed == 1);
		}
		// implementation of ran_likelihood
		template <typename Vector>
		Vector template_ran_likelihood(
			const Vector& theta  ,
			const Vector& u      )
		{	typedef typename Vector::value_type scalar;

			Vector vec(1);

			// initialize part of log-density that is always smooth
			vec[0] = scalar(0.0);

			// p(y | u, theta) , p(u_i | theta)
			for(size_t i = 0; i < y_.size(); ++i)
			{	scalar res = y_[i] - theta[0] * exp(u[i]);
				vec[0] += res * res / scalar(2.0);
				vec[0] += u[i] * u[i] / scalar(2.0);
			}

			return vec;
		}
		// a1_vector version of ran_likelihood
		virtual a1_vector ran_likelihood(
			const a1_vector& fixed_vec, const a1_vector& random_vec
		)
		{	return template_ran_likelihood( fixed_vec, random_vec ); }
	};
}

bool laplace_obj_tst(void)
{
	bool   ok = true;

	size_t n_data   = 1;
	size_t n_fixed  = 1;
	size_t n_random = n_data;
	d_vector y(n_data), fixed_in(n_fixed), random_in(n_random);
	d_vector fixed_lower(n_fixed), fixed_upper(n_fixed);
	d_vector random_lower(n_random), random_upper(n_random);

	double inf = std::numeric_limits<double>::infinity();
	fixed_in[0]     = 1.0;
	fixed_lower[0]  = 1e-5;
	fixed_upper[0]  = 1e+2;
	for(size_t i = 0; i < n_data; ++i)
	{	random_lower[i] = -inf;
		random_upper[i] = +inf;
		random_in[i]    = 0.0;
		y[i]            = double(20 + i);
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed   = false;
	bool bool_sparsity = true;
	CppAD::mixed::d_sparse_rcv A_rcv; // empty matrix
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, y
	);

	// initialize mixed_object
	mixed_object.initialize(fixed_in, random_in);

	// optimize with respect to fixed effects

	// optimize the fixed effects using quasi-Newton method
	std::string fixed_ipopt_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           first-order\n"
		"String  derivative_test_print_all yes\n"
		"Numeric tol                       1e-10\n"
		"Integer max_iter                  30\n"
	;
	std::string random_ipopt_options =
		"Integer print_level 0\n"
		"String  sb          yes\n"
		"String  derivative_test second-order\n"
	;
	vector<double> fixed_scale = fixed_in;
	vector<double> fix_constraint_lower(0), fix_constraint_upper(0);
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
	vector<double> fixed_opt = solution.fixed_opt;

	// optimize the random effects
	vector<double> random_opt = mixed_object.optimize_random(
		random_ipopt_options, fixed_opt, random_lower, random_upper, random_in
	);

	// check that derivative is zero
	vector<double> r_fixed(n_fixed);
	mixed_object.ran_obj_jac(fixed_opt, random_opt, r_fixed);

	//
	for(size_t i = 0; i < n_fixed; ++i)
		ok &= std::fabs(r_fixed[i]) < 1e-6;
	//

	return ok;
}
// END C++
