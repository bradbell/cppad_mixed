// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
// test of using ran_likelihood_hes for computation of random Hessian
/*
$head Model$$
$latex \[
	\B{p}( y_i | \theta , u ) \sim \B{N} ( \exp(u_i) + \theta_0 , \theta_1^2 )
\] $$
$latex \[
	\B{p}( u_i | \theta ) \sim \B{N} ( 0 , 1 )
\] $$
$latex \[
	\B{p}( \theta ) \sim \B{N} ( 4 , 1 )
\] $$
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>

namespace {
	using CppAD::vector;
	using CppAD::log;
	using CppAD::AD;
	using CppAD::mixed::sparse_mat_info;

	// use the default ran_likelihood_hes
	bool default_ran_likelihood_hes_ = false;

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
		// ------------------------------------------------------------------
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

			// pi
			Float sqrt_2pi = Float( CppAD::sqrt(8.0 * CppAD::atan(1.0) ) );

			for(size_t i = 0; i < n_random_; i++)
			{	Float mu     = exp( u[i] ) * theta[0];
				Float sigma  = theta[1];
				Float res    = (y_[i] - mu) / sigma;

				// p(y_i | u, theta)
				vec[0] += log(sqrt_2pi * sigma) + res * res / Float(2.0);

				// p(u_i | theta)
				vec[0] += log(sqrt_2pi) + u[i] * u[i] / Float(2.0);
			}
			return vec;
		}
		// ------------------------------------------------------------------
		// ran_likelihood_hes
		vector<a1_double> ran_likelihood_hes(
			const vector<a1_double>& theta  ,
			const vector<a1_double>& u      ,
			const vector<size_t>&    row    ,
			const vector<size_t>&    col    )
		{	//
			if( default_ran_likelihood_hes_ )
				return vector<a1_double>(0);
			//
			assert( col.size() == row.size() );

			// return value
			vector<a1_double> val(n_random_);

			// for each component of the return value
			for(size_t k = 0; k < n_random_; k++)
			{	// initialize it as zero
				val[k] = a1_double(0.0);

				// for this ran_likelihood only the diagonal is non-zero
				if( row[k] == col[k] )
				{	size_t i = row[k];
					//
					a1_double mu        = exp( u[i] ) * theta[0];
					a1_double sigma     = theta[1];
					a1_double res       = (y_[i] - mu) / sigma;
					a1_double res_ui    = - mu / sigma;
					a1_double res_ui_ui = - mu / sigma;
					a1_double sq_ui     = res * res_ui;
					a1_double sq_ui_ui  = res_ui * res_ui + res * res_ui_ui;
					val[k]              = sq_ui_ui + 1.0;
				}
			}
			return val;
		}
		// ------------------------------------------------------------------
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
}

bool ran_likelihood_hes(void)
{
	bool   ok = true;
	double inf = std::numeric_limits<double>::infinity();
	double tol = 1e-8;

	size_t n_data   = 10;
	size_t n_fixed  = 2;
	size_t n_random = n_data;
	vector<double>
		fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
	fixed_lower[0] = - inf; fixed_in[0] = 5.0; fixed_upper[0] = inf;
	fixed_lower[1] = .01;   fixed_in[1] = 0.5; fixed_upper[1] = inf;
	//
	// explicit constriants (in addition to l1 terms)
	vector<double> fix_constraint_lower(0), fix_constraint_upper(0);
	//
	vector<double> data(n_data), random_in(n_random);
	for(size_t i = 0; i < n_data; i++)
	{	data[i]       = double(i + 1);
		random_in[i] = 1.0;
	}

	// Newton method uses more derivatives than quasi-newton method
	bool quasi_fixed = false;

	// random constraint matrix
	CppAD::mixed::sparse_mat_info A_info; // empty matrix

	// optimize the fixed effects using quasi-Newton method
	// If the derivatives are correct, the optimzation converges in 6
	// iterations. If convergence fails, change print_level to 5
	std::string fixed_ipopt_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           second-order\n"
		"String  derivative_test_print_all no\n"
		"Numeric tol                       1e-8\n"
		"Integer max_iter                  7\n"
	;
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
	vector< vector<double> > theta_out(2);
	for(size_t user_defined = 0; user_defined < 2; user_defined++)
	{	default_ran_likelihood_hes_ = ! default_ran_likelihood_hes_;

		// object that is derived from cppad_mixed
		mixed_derived mixed_object(
			n_fixed, n_random, quasi_fixed, A_info, data
		);
		mixed_object.initialize(fixed_in, random_in);
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
		theta_out[user_defined] = solution.fixed_opt;
	}

	// check that the results were the same
	// with and without the user defined ran_likelihood_hes
	double eps = 5.0 * tol;
	ok &= CppAD::NearEqual( theta_out[0][0], theta_out[1][0], eps, eps );
	ok &= CppAD::NearEqual( theta_out[0][1], theta_out[1][1], eps, eps );

	return ok;
}
// END C++
