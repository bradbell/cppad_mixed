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
$begin opt_ran_nan.cpp$$
$spell
	CppAD
	cppad
$$

$section Nan's During Optimization of Random Effects: Example and Test$$

$code
$srcfile%example/user/opt_ran_nan.cpp
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
	using CppAD::mixed::a2_double;
	using CppAD::mixed::d_vector;
	using CppAD::mixed::a2_vector;
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
			const sparse_rcv&      A_rcv         ,
			const d_vector&       y              ) :
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

			// initialize part of log-density that is always smooth
			vec[0] = scalar(0.0);

			// sqrt_2pi = CppAD::sqrt( 8.0 * CppAD::atan(1.0) );

			// sum of residual squared
			scalar sum_sq  = scalar(0.0);
			for(size_t i = 0; i < y_.size(); i++)
			{	scalar mu     = u[i];
				scalar sigma  = theta[i];
				scalar res    = (y_[i] - mu) / sigma;

				// Gaussian likelihood
				vec[0]  += log(sigma) + res * res / scalar(2.0);
				// following term does not depend on fixed or random effects
				// vec[0]  += log(sqrt_2pi);

				// add to sum
				sum_sq += res * res;
			}

			// return nan when sum of squares is less than 1e-4
			scalar vec_0 = CppAD::numeric_limits<scalar>::quiet_NaN();
			scalar small = scalar( 1e-4 );
			vec_0  = CppAD::CondExpGt(sum_sq, + small, vec[0], vec_0);
			vec[0] = vec_0;
			//
			return vec;
		}
		// a2_vector version of ran_likelihood
		virtual a2_vector ran_likelihood(
			const a2_vector& fixed_vec, const a2_vector& random_vec
		)
		{	return template_ran_likelihood( fixed_vec, random_vec ); }
		// a3_vector version of ran_likelihood
		virtual a3_vector ran_likelihood(
			const a3_vector& fixed_vec, const a3_vector& random_vec
		)
		{	return template_ran_likelihood( fixed_vec, random_vec ); }
		// we expect to get a warnings
		virtual void warning(const std::string& warning_message)
		{ }
	};
}

bool opt_ran_nan_xam(void)
{
	bool   ok = true;

	size_t n_data = 10;
	d_vector data(n_data), fixed_vec(n_data), random_in(n_data);

	for(size_t i = 0; i < n_data; i++)
	{	data[i]      = double(i + 1);
		fixed_vec[i] = 1.0;
		random_in[i] = 0.0;
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed   = true;
	bool bool_sparsity = true;
	sparse_rcv A_rcv; // empty matrix
	mixed_derived mixed_object(
		n_data, n_data, quasi_fixed, bool_sparsity, A_rcv, data
	);
	mixed_object.initialize(fixed_vec, random_in);

	// lower and upper limits for random effects
	double inf = std::numeric_limits<double>::infinity();
	d_vector random_lower(n_data), random_upper(n_data);
	for(size_t i = 0; i < n_data; i++)
	{	random_lower[i] = -inf;
		random_upper[i] = +inf;
	}

	// -----------------------------------------------------------------------
	// attempt to use ipopt to determine the optimal random effects
	std::string ipopt_options;
	ipopt_options += "Integer print_level      0\n";
	ipopt_options += "Integer max_iter         10\n";
	ipopt_options += "String  sb               yes\n";
	ipopt_options += "String  derivative_test  second-order\n";
	d_vector random_out = mixed_object.optimize_random(
		ipopt_options, fixed_vec, random_lower, random_upper, random_in
	);

	// check that the optimize had backed up and solve problem
	// to be near forbidden region (where nans occur)
	double sum_sq = 0.0;
	for(size_t i = 0; i < n_data; i++)
	{	double mu    = random_out[i];
		double sigma = fixed_vec[i];
		double res   = (data[i] - mu) / sigma;
		sum_sq      += res * res;
	}
	ok &= sum_sq >= 1e-4;
	ok &= sum_sq <= 1e-3;
	//
	return ok;
}
// END C++
