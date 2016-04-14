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
$section Sample From Fixed Effects Posterior: Example and Test$$
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/sparse_mat_info.hpp>
# include <cppad/mixed/manage_gsl_rng.hpp>
# include <Eigen/Dense>

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

			// pi
			Float sqrt_2pi = Float( CppAD::sqrt(8.0 * CppAD::atan(1.0) ) );

			for(size_t i = 0; i < n_random_; i++)
			{	Float mu     = u[i] + theta[0];
				Float sigma  = theta[1];
				Float res    = (y_[i] - mu) / sigma;

				// p(y_i | u, theta)
				vec[0] += log(sqrt_2pi * sigma) + res * res / Float(2.0);

				// p(u_i | theta)
				vec[0] += log(sqrt_2pi) + u[i] * u[i] / Float(2.0);
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
}

bool sample_fixed(void)
{
	bool   ok = true;
	double inf = std::numeric_limits<double>::infinity();
	//
	// initialize gsl random number generator
	size_t random_seed = CppAD::mixed::new_gsl_rng(0);
	//
	size_t n_data   = 10;
	size_t n_fixed  = 2;
	size_t n_random = n_data;
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
	bool quasi_fixed = true;
	CppAD::mixed::sparse_mat_info A_info; // empty matrix
	mixed_derived mixed_object(n_fixed, n_random, quasi_fixed, A_info, data);
	mixed_object.initialize(fixed_in, random_in);

	// optimize the fixed effects using quasi-Newton method
	std::string fixed_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           first-order\n"
		"String  derivative_test_print_all yes\n"
		"Numeric tol                       1e-8\n"
	;
	std::string random_options =
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
	// optimize fixed effects
	CppAD::mixed::fixed_solution solution = mixed_object.optimize_fixed(
		fixed_options,
		random_options,
		fixed_lower,
		fixed_upper,
		fix_constraint_lower,
		fix_constraint_upper,
		fixed_in,
		random_lower,
		random_upper,
		random_in
	);
	//
	// now change the limit on the second fixed effect so it is active
	fixed_lower[1] = 1.5 * solution.fixed_opt[1];
	fixed_in[1]    = 2.0 * solution.fixed_opt[1];
	//
	// optimize fixed effects
	solution = mixed_object.optimize_fixed(
		fixed_options,
		random_options,
		fixed_lower,
		fixed_upper,
		fix_constraint_lower,
		fix_constraint_upper,
		fixed_in,
		random_lower,
		random_upper,
		random_in
	);
	//
	// check that one of the constraints is active
	ok &= solution.fixed_lag.size() == n_fixed;
	ok &= solution.fixed_lag[0] == 0.0;
	ok &= solution.fixed_lag[1] != 0.0;
	ok &= solution.fix_con_lag.size() == 0;
	ok &= solution.ran_con_lag.size() == 0.0;
	//
	// corresponding optimal random effects
	vector<double> random_opt = mixed_object.optimize_random(
		random_options,
		solution.fixed_opt,
		random_lower,
		random_upper,
		random_in
	);
	//
	// compute corresponding information matrix
	CppAD::mixed::sparse_mat_info
	information_info = mixed_object.information_mat(solution, random_opt);
	//
	// sample from the posterior for fixed effects
	size_t n_sample = 10000;
	CppAD::vector<double> sample( n_sample * n_fixed );
	mixed_object.sample_fixed(
		sample,
		information_info,
		solution,
		random_opt
	);
	//
	typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > matrix;
	//
	// compute sample covariance matrix
	matrix sample_cov = matrix::Zero(n_fixed, n_fixed);
	for(size_t i = 0; i < n_sample; i++)
	{	matrix diff(n_fixed, 1);
		for(size_t j = 0; j < n_fixed; j++)
			diff(j, 0) = sample[ i * n_fixed + j] - solution.fixed_opt[j];
		sample_cov += diff * diff.transpose();
	}
	sample_cov *= 1.0 / double(n_sample);
	//
	matrix info_mat(n_fixed, n_fixed);
	size_t K = ( n_fixed * (n_fixed + 1) ) / 2;
	ok &= K == information_info.row.size();
	for(size_t k = 0; k < K; k++)
	{	size_t i = information_info.row[k];
		size_t j = information_info.col[k];
		info_mat(i, j) = information_info.val[k];
		info_mat(j, i) = information_info.val[k];
	}
	matrix cov_mat = info_mat.inverse();
	//
	for(size_t i = 0; i < n_fixed; i++)
	{	for(size_t j = 0; j < n_fixed; j++)
		{	double value = sample_cov(i, j);
			double check = 0.0;
			if( i == 0 && j == 0 )
				check = cov_mat(i, j);
			double scale = std::sqrt( cov_mat(i, i) * cov_mat(j, j) );
			ok &= std::fabs(value - check) / scale < .05;
		}
	}
	//
	if( ! ok )
		std::cout << "\nrandom_seed = " << random_seed << "\n";
	//
	return ok;
}
// END C++
