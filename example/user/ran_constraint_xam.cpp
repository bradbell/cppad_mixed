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
$begin ran_constraint_xam.cpp$$
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
$srcfile%example/user/ran_constraint_xam.cpp
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

	class mixed_derived : public cppad_mixed {
	private:
		const size_t          n_fixed_;
		const size_t          n_random_;
		const vector<double>& y_;
	// ----------------------------------------------------------------------
	public:
		// constructor
		mixed_derived(
			size_t n_fixed                    ,
			size_t n_random                   ,
			bool   quasi_fixed                ,
			const  sparse_mat_info& A_info    ,
			const vector<double>& y           ) :
			cppad_mixed(n_fixed, n_random, quasi_fixed, A_info) ,
			n_fixed_(n_fixed)     ,
			n_random_(n_random)   ,
			y_(y)
		{	assert( n_fixed == 2);
			assert( y_.size() == n_random_ );
		}
	// ----------------------------------------------------------------------
	private:
		// implementation of ran_likelihood
		template <class Float>
		vector<Float> implement_ran_likelihood(
			const vector<Float>& theta  ,
			const vector<Float>& u      )
		{	assert( theta.size() == n_fixed_ );
			assert( u.size() == y_.size() );
			vector<Float> vec(1);

			// initialize part of log-density that is always smooth
			vec[0] = Float(0.0);

			// Float sqrt_2pi = Float( CppAD::sqrt(8.0 * CppAD::atan(1.0) ) );

			for(size_t i = 0; i < n_random_; i++)
			{	Float mu     = u[i] + theta[0];
				Float sigma  = theta[1];
				Float res    = (y_[i] - mu) / sigma;

				// p(y_i | u, theta)
				vec[0] += log(sigma) + res * res / Float(2.0);
				// following term does not depend on fixed or random effects
				// vec[0] += log(sqrt_2pi);

				// p(u_i | theta)
				vec[0] += u[i] * u[i] / Float(2.0);
				// following term does not depend on fixed or random effects
				// vec[0] += log(sqrt_2pi);
			}
			return vec;
		}
		// implementation of fix_likelihood
		template <class Float>
		vector<Float> implement_fix_likelihood(
			const vector<Float>& fixed_vec  )
		{	assert( fixed_vec.size() == n_fixed_ );
			vector<Float> vec(1);

			// initialize part of log-density that is smooth
			vec[0] = Float(0.0);

			// compute these factors once
			Float sqrt_2pi = Float( CppAD::sqrt( 8.0 * CppAD::atan(1.0) ) );

			for(size_t j = 0; j < n_fixed_; j++)
			{	Float mu     = Float(4.0);
				Float sigma  = Float(1.0);
				Float res    = (fixed_vec[j] - mu) / sigma;

				// This is a Gaussian term, so entire density is smooth
				vec[0]  += log(sqrt_2pi * sigma) + res * res / Float(2.0);
			}
			return vec;
		}
	// ----------------------------------------------------------------------
	public:
		// User defined virtual functions
		virtual vector<a2_double> ran_likelihood(
			const vector<a2_double>& fixed_vec  ,
			const vector<a2_double>& random_vec )
		{	return implement_ran_likelihood(fixed_vec, random_vec); }
		virtual vector<a1_double> ran_likelihood(
			const vector<a1_double>& fixed_vec  ,
			const vector<a1_double>& random_vec )
		{	return implement_ran_likelihood(fixed_vec, random_vec); }
		//
		virtual vector<a1_double> fix_likelihood(
			const vector<a1_double>& fixed_vec  )
		{	return implement_fix_likelihood(fixed_vec); }
		// ------------------------------------------------------------------
	};
	double sum_random_effects(
		size_t n_random, const CppAD::mixed::sparse_mat_info& A_info
	)
	{
		double inf = std::numeric_limits<double>::infinity();

		size_t n_fixed  = 2;
		size_t n_data   = n_random;
		vector<double>
			fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
		fixed_lower[0] = - inf; fixed_in[0] = 2.0; fixed_upper[0] = inf;
		fixed_lower[1] = .01;   fixed_in[1] = 0.5; fixed_upper[1] = inf;
		//
		// explicit constriants (in addition to l1 terms)
		vector<double> fix_constraint_lower(0), fix_constraint_upper(0);
		//
		vector<double> data(n_data), random_in(n_random);
		for(size_t i = 0; i < n_data; i++)
		{	data[i]       = double(i + 1);
			random_in[i] = 0.0;
		}

		// object that is derived from cppad_mixed
		bool quasi_fixed = false;
		mixed_derived mixed_object(
			n_fixed, n_random, quasi_fixed, A_info, data
		);
		mixed_object.initialize(fixed_in, random_in);

		// optimize the fixed effects using quasi-Newton method
		std::string fixed_ipopt_options =
			"Integer print_level               0\n"
			"String  sb                        yes\n"
			"String  derivative_test           second-order\n"
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
		vector<double> random_lower(n_random), random_upper(n_random);
		for(size_t i = 0; i < n_random; i++)
		{	random_lower[i] = -inf;
			random_upper[i] = +inf;
		}
		// optmize fixed effects
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
		// corresponding optimal random effects
		vector<double> random_out = mixed_object.optimize_random(
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
	CppAD::mixed::sparse_mat_info A_info;
	double sum = sum_random_effects(n_random, A_info);
	ok        &= CppAD::abs(sum) > 0.5;

	// constrain sum of random effects to be zero
	A_info.resize(n_random);
	for(size_t k = 0; k < n_random; k++)
	{	A_info.row[k] = 0;
		A_info.col[k] = k;
		A_info.val[k] = 1.0;
	}
	sum = sum_random_effects(n_random, A_info);
	ok &= CppAD::abs(sum) < 2.0 * tol;

	return ok;
}
// END C++
