// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-17 University of Washington
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
	//
	using CppAD::mixed::sparse_rcv;
	using CppAD::mixed::a1_double;
	using CppAD::mixed::a2_double;
	using CppAD::mixed::d_vector;
	using CppAD::mixed::a1_vector;
	using CppAD::mixed::a2_vector;
	//
	class mixed_derived : public cppad_mixed {
	private:
		const d_vector&       y_;
		bool                  ran_likelihood_jac_called_;
	public:
		// constructor
		mixed_derived(
			size_t                 n_fixed       ,
			size_t                 n_random      ,
			bool                   quasi_fixed   ,
			bool                   bool_sparsity ,
			const sparse_rcv&      A_rcv         ,
			const d_vector&        y             ) :
			cppad_mixed(
				n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
			),
			y_(y),
			ran_likelihood_jac_called_(false)
		{ }
		// ------------------------------------------------------------------
		// implementation of ran_likelihood
		template <typename Vector>
		Vector template_ran_likelihood(
			const Vector&         theta  ,
			const Vector&         u      )
		{	typedef typename Vector::value_type scalar;

			Vector vec(1);

			// sqrt_2pi = CppAD::sqrt( 8.0 * CppAD::atan(1.0) );

			// initialize summation
			vec[0] = scalar(0.0);

			// for each data and random effect
			for(size_t i = 0; i < y_.size(); i++)
			{	// data residual
				scalar model = theta[i] + u[i];
				scalar res = (y_[i] - model);

				// This is a Gaussian term, so density is smooth
				vec[0]  += res * res / scalar(2.0);

				// prior for random effects
				vec[0]  += u[i] * u[i] / scalar(2.0);
			}
			return vec;
		}
		// a2_vector version of ran_likelihood
		// a3_vector version of ran_likelihood
		virtual a3_vector ran_likelihood(
			const a3_vector& fixed_vec, const a3_vector& random_vec
		)
		{	return template_ran_likelihood( fixed_vec, random_vec ); }
		// ------------------------------------------------------------------
		// ran_likelihood_jac
		template <typename Vector>
		Vector template_ran_likelihood_jac(
			const Vector&         theta  ,
			const Vector&         u      )
		{	typedef typename Vector::value_type scalar;

			// return value
			Vector vec(y_.size());

			// for each data and random effect
			for(size_t i = 0; i < y_.size(); i++)
			{	scalar model  = theta[i] + u[i];
				scalar res    = y_[i] - model;
				scalar res_ui;
				res_ui        = - 1.0 / 2.0;

				// This is a Gaussian term, so entire density is smooth
				vec[i]  =  2.0 * res * res_ui;
				vec[i] +=  u[i];
			}
			return vec;
		}
		// a1_vector version of ran_likelihood_jac
		virtual a2_vector ran_likelihood_jac(
			const a2_vector& theta, const a2_vector& u
		)
		{	ran_likelihood_jac_called_ = true;
			return template_ran_likelihood_jac( theta, u );
		}
		bool ran_likelihood_jac_called(void) const
		{	return ran_likelihood_jac_called_; }
	};
}

bool ran_likelihood_jac_xam(void)
{
	bool   ok  = true;
	double inf = std::numeric_limits<double>::infinity();
	//
	size_t n_data   = 10;
	size_t n_fixed  = n_data;
	size_t n_random = n_data;
	d_vector data(n_data);
	d_vector fixed_lower(n_data), fixed_in(n_data), fixed_upper(n_data);
	d_vector random_lower(n_data), random_in(n_data), random_upper(n_data);
	for(size_t i = 0; i < n_data; i++)
	{	data[i]       = double(i + 1);
		//
		fixed_lower[i] = random_lower[i] = -inf;
		fixed_upper[i] = random_upper[i] = +inf;
		fixed_in[i]    = random_in[i]    = 1.0;
	}
	//
	// object that is derived from cppad_mixed
	bool quasi_fixed   = true;
	bool bool_sparsity = true;
	sparse_rcv A_rcv; // empty matrix
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
		"Numeric tol             1e-10\n"
	;
	//
	for(size_t i = 0; i < n_random; i++)
	{	random_lower[i] = -inf;
		random_upper[i] = +inf;
	}
	// ------------------------------------------------------------------
	// optimize fixed effects
	d_vector fixed_scale = fixed_in;
	d_vector fix_constraint_lower(0), fix_constraint_upper(0);
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

	for(size_t i = 0; i < n_data; i++)
		ok &= CppAD::NearEqual( data[i], fixed_out[i], 1e-7, 1e-7);
	ok &= mixed_object.ran_likelihood_jac_called();


	return ok;
}
// END C++
