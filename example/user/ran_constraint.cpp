// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin ran_constraint.cpp$$
$spell
	CppAD
	cppad
	hes
	eval
	interp
	xam
$$

$section Constraints On Random Effects: Example and Test$$

This example demonstrates
$cref/random constraints/cppad_mixed/Problem/Random Constraints/$$.
To be specific, it demonstrates a case where the constraints ensure
that the sum of the
$cref/optional random effects
	/cppad_mixed
	/Notation
	/Optimal Random Effects, u^(theta)
/$$
is zero.
In addition, for the same case without the constraint,
the optimal random effects do not satisfy this condition.

$code
$srcfile%example/user/ran_constraint.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>

namespace {
	using CppAD::log;
	using CppAD::AD;
	//
	using CppAD::mixed::d_sparse_rcv;
	using CppAD::mixed::a1_double;
	using CppAD::mixed::d_vector;
	using CppAD::mixed::a1_vector;
	//
	class mixed_derived : public cppad_mixed {
	private:
		const size_t          n_fixed_;
		const size_t          n_random_;
		const d_vector&       y_;
	// ----------------------------------------------------------------------
	public:
		// constructor
		mixed_derived(
			size_t                 n_fixed       ,
			size_t                 n_random      ,
			bool                   quasi_fixed   ,
			bool                   bool_sparsity ,
			const d_sparse_rcv&    A_rcv         ,
			const d_vector&       y              ) :
			cppad_mixed(
				n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
			)                     ,
			n_fixed_(n_fixed)     ,
			n_random_(n_random)   ,
			y_(y)
		{	assert( n_fixed == 2);
			assert( y_.size() == n_random_ );
		}
		// implementation of ran_likelihood
		template <typename Vector>
		Vector template_ran_likelihood(
			const Vector&         theta  ,
			const Vector&         u      )
		{	typedef typename Vector::value_type scalar;

			assert( theta.size() == n_fixed_ );
			assert( u.size() == y_.size() );
			Vector vec(1);

			// initialize part of log-density that is always smooth
			vec[0] = scalar(0.0);

			// sqrt_2pi = CppAD::sqrt(8.0 * CppAD::atan(1.0) );

			for(size_t i = 0; i < n_random_; i++)
			{	scalar mu     = u[i] + theta[0];
				scalar sigma  = theta[1];
				scalar res    = (y_[i] - mu) / sigma;

				// p(y_i | u, theta)
				vec[0] += log(sigma) + res * res / scalar(2.0);
				// following term does not depend on fixed or random effects
				// vec[0] += log(sqrt_2pi);

				// p(u_i | theta)
				vec[0] += u[i] * u[i] / scalar(2.0);
				// following term does not depend on fixed or random effects
				// vec[0] += log(sqrt_2pi);
			}
			return vec;
		}
		// a1_vector version of ran_likelihood
		virtual a1_vector ran_likelihood(
			const a1_vector& fixed_vec, const a1_vector& random_vec
		)
		{	return template_ran_likelihood( fixed_vec, random_vec ); }
		// implementation of fix_likelihood
		template <typename Vector>
		Vector template_fix_likelihood(
			const Vector&         fixed_vec  )
		{	typedef typename Vector::value_type scalar;

			assert( fixed_vec.size() == n_fixed_ );
			Vector vec(1);

			// initialize part of log-density that is smooth
			vec[0] = scalar(0.0);

			// compute these factors once
			scalar sqrt_2pi = scalar(
				 CppAD::sqrt( 8.0 * CppAD::atan(1.0)
			));

			for(size_t j = 0; j < n_fixed_; j++)
			{	scalar mu     = scalar(4.0);
				scalar sigma  = scalar(1.0);
				scalar res    = (fixed_vec[j] - mu) / sigma;

				// This is a Gaussian term, so entire density is smooth
				vec[0]  += log(sqrt_2pi * sigma) + res * res / scalar(2.0);
			}
			return vec;
		}
		// a1_vector version of fix_likelihood
		virtual a1_vector fix_likelihood(const a1_vector& fixed_vec)
		{	return template_fix_likelihood( fixed_vec ); }
	};
	// ----------------------------------------------------------------------
	double sum_random_effects(
		size_t n_random, const d_sparse_rcv&    A_rcv
	)
	{
		double inf = std::numeric_limits<double>::infinity();

		size_t n_fixed  = 2;
		size_t n_data   = n_random;
		d_vector
			fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
		fixed_lower[0] = - inf; fixed_in[0] = 2.0; fixed_upper[0] = inf;
		fixed_lower[1] = .01;   fixed_in[1] = 0.5; fixed_upper[1] = inf;
		//
		// explicit constriants (in addition to l1 terms)
		d_vector fix_constraint_lower(0), fix_constraint_upper(0);
		//
		d_vector data(n_data), random_in(n_random);
		for(size_t i = 0; i < n_data; i++)
		{	data[i]       = double(i + 1);
			random_in[i] = 0.0;
		}

		// object that is derived from cppad_mixed
		bool quasi_fixed   = false;
		bool bool_sparsity = false;
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
			"Numeric tol                       1e-8\n"
		;
		// random_ipopt_options is non-empty, so using ipopt for random effects
		std::string random_ipopt_options =
			"Integer print_level     0\n"
			"String  sb              yes\n"
			"String  derivative_test second-order\n"
			"Numeric tol             1e-8\n"
		;
		d_vector random_lower(n_random), random_upper(n_random);
		for(size_t i = 0; i < n_random; i++)
		{	random_lower[i] = -inf;
			random_upper[i] = +inf;
		}
		// optmize fixed effects
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
		//
		// corresponding optimal random effects
		d_vector random_out = mixed_object.optimize_random(
			random_ipopt_options,
			fixed_out,
			random_lower,
			random_upper,
			random_in
		);

		// compute return value
		double sum = 0.0;
		for(size_t i = 0; i < n_random; i++)
			sum += random_out[i];
		//
		return sum;
	}
}
bool ran_constraint_xam(void)
{	bool ok         = true;
	double tol      = 1e-8;
	size_t n_random = 10;

	// empty matrix (no constraints)
	d_sparse_rcv A_empty;
	double sum = sum_random_effects(n_random, A_empty);
	ok        &= fabs(sum) > 0.5;

	// constrain sum of random effects to be zero
	CppAD::mixed::sparse_rc A_pattern(1, n_random, n_random);
	for(size_t k = 0; k < n_random; k++)
		A_pattern.set(k, 0, k);
	d_sparse_rcv A_rcv(A_pattern);
	for(size_t k = 0; k < n_random; k++)
		A_rcv.set(k, 1.0);
	sum = sum_random_effects(n_random, A_rcv);
	ok &= fabs(sum) < 2.0 * tol;

	return ok;
}
// END C++
