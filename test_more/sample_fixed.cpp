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
	typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > double_mat;
	//
	// Used for debugging
	// void print(const char* name , const double_mat& mat)
	// {	std::cout << "\n" << name << " =\n" << mat << "\n"; }
	//
	class mixed_derived : public cppad_mixed {
	private:
		size_t                n_fixed_;
		const vector<double>& y_;
	public:
		// constructor
		mixed_derived(
			size_t n_fixed                    ,
			size_t n_random                   ,
			bool   quasi_fixed                ,
			const  sparse_mat_info& A_info    ,
			const vector<double>& y           ) :
			cppad_mixed(n_fixed, n_random, quasi_fixed, A_info) ,
			n_fixed_(n_fixed)                                      ,
			y_(y)
		{}
	private:
		// implementation of ran_likelihood
		template <class Float>
		vector<Float> implement_ran_likelihood(
			const vector<Float>& theta  ,
			const vector<Float>& u      )
		{	assert( u.size() == y_.size() );
			assert( theta.size() == y_.size() );
			vector<Float> vec(1);

			// initialize part of log-density that is always smooth
			vec[0] = Float(0.0);

			// pi
			Float sqrt_2pi = Float( CppAD::sqrt(8.0 * CppAD::atan(1.0) ) );

			for(size_t i = 0; i < y_.size(); i++)
			{	Float mu     = u[i] + theta[i] + theta[0];
				Float sigma  = Float(1.0);
				Float res    = (y_[i] - mu) / sigma;

				// p(y_i | u, theta)
				vec[0] += log(sqrt_2pi * sigma) + res*res / Float(2.0);

				// p(u_i | theta)
				vec[0] += log(sqrt_2pi) + u[i] * u[i] / Float(2.0);
			}
			return vec;
		}
		// ------------------------------------------------------------------
		// fix_constraint
		template <class Float>
		vector<Float> implement_fix_constraint(
			const vector<Float>& fixed_vec  )
		{	vector<Float> ret_val(2);
			//
			for(size_t i = 0; i < 2; i++)
			{	ret_val[i] = 0.0;
				for(size_t j = 0; j < fixed_vec.size(); j++)
				{	if( i == 0 )
						ret_val[i] += fixed_vec[j];
					else
						ret_val[i] += fixed_vec[j] * fixed_vec[j];
				}
			}
			return ret_val;
		}
	public:
		// ------------------------------------------------------------------
		// ran_likelihood
		virtual vector<a2_double> ran_likelihood(
			const vector<a2_double>& fixed_vec  ,
			const vector<a2_double>& random_vec )
		{	return implement_ran_likelihood(fixed_vec, random_vec); }
		virtual vector<a1_double> ran_likelihood(
			const vector<a1_double>& fixed_vec  ,
			const vector<a1_double>& random_vec )
		{	return implement_ran_likelihood(fixed_vec, random_vec); }
		// ------------------------------------------------------------------
		// fix_constraint
		virtual vector<a1_double> fix_constraint(
			const vector<a1_double>& fixed_vec  )
		{	return implement_fix_constraint(fixed_vec); }
		vector<double> fix_constraint(
			const vector<double>& fixed_vec  )
		{	return implement_fix_constraint(fixed_vec); }
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
	size_t n_data   = 4;
	size_t n_fixed  = n_data;
	size_t n_random = n_data;
	vector<double>
		fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
	for(size_t i = 0; i < n_fixed; i++)
	{	fixed_lower[i] = - inf;
		fixed_in[i]    = 0.1;
		fixed_upper[i] = inf;
	}
	//
	// explicit constriants (in addition to l1 terms)
	vector<double> fix_constraint_lower(2), fix_constraint_upper(2);
	for(size_t i = 0; i < 2; i++)
	{	fix_constraint_lower[i] = -inf;
		fix_constraint_upper[i] = inf;
	}
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

	// options for optimization
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
	// lower and uppser limits for random effects
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
	// now change the limit on the first fixed effect so it is active
	size_t bnd_j   = 0;
	fixed_lower[bnd_j] = 1.5 * solution.fixed_opt[bnd_j];
	fixed_in[bnd_j]    = 2.0 * solution.fixed_opt[bnd_j];
	//
	// also change the limit on the first fix_constraint so it is active
	vector<double> fix_con = mixed_object.fix_constraint(solution.fixed_opt);
	size_t con_i = 0;
	fix_constraint_lower[con_i] = 1.5 * fix_con[con_i];
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
	for(size_t j = 0; j < n_fixed; j++)
	{	if( j == bnd_j )
			ok &= solution.fixed_lag[j] != 0.0;
		else
			ok &= solution.fixed_lag[j] == 0.0;
	}
	ok &= solution.fix_con_lag.size() == 2;
	for(size_t i = 0; i < 2; i++)
	{	if( i == con_i )
			ok &= solution.fix_con_lag[i] != 0.0;
		else
			ok &= solution.fix_con_lag[i] == 0.0;
	}
	ok &= solution.ran_con_lag.size() == 0;
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
	// compute sample covariance matrix
	double_mat sample_cov = double_mat::Zero(n_fixed, n_fixed);
	for(size_t i = 0; i < n_sample; i++)
	{	double_mat diff(n_fixed, 1);
		for(size_t j = 0; j < n_fixed; j++)
			diff(j, 0) = sample[ i * n_fixed + j] - solution.fixed_opt[j];
		sample_cov += diff * diff.transpose();
	}
	sample_cov *= 1.0 / double(n_sample);
	//
	// unconstrained information matrix
	double_mat full_info = double_mat::Zero(n_fixed, n_fixed);
	size_t K = information_info.row.size();
	for(size_t k = 0; k < K; k++)
	{	size_t i = information_info.row[k];
		size_t j = information_info.col[k];
		full_info(i, j) = information_info.val[k];
		full_info(j, i) = information_info.val[k];
	}
	//
	// consider the first and second fixed effects to be functions
	// of the other fixed effect as follows
	// theta[0] = lower[0]
	// theta[1] = fix_con[0] - lower[0] - theta[2] - ...
	//
	// theta[0:end] = E * theta[2:end] + b
	//
	double_mat E = double_mat::Zero(n_fixed, n_fixed - 2);
	for(size_t j = 0; j < n_fixed - 2; j++)
		E(1, j) = - 1.0;
	for(size_t i = 2; i < n_fixed; i++)
		E(i, i-2) = 1.0;
	//
	// information matrix on the reduced set of variables
	double_mat reduced_info = E.transpose() * full_info * E;
	//
	// covariance on the reduced set of variables
	double_mat reduced_cov  = reduced_info.inverse();
	//
	// covariance on the full set of variables
	double_mat full_cov = E * reduced_cov * E.transpose();
	//
	double max_cov = 0.0;
	for(size_t i = 0; i < n_fixed; i++)
		max_cov = std::max( max_cov, full_cov(i, i) );
	//
	//
	double max_err = 0.0;
	for(size_t i = 0; i < n_fixed; i++)
	{	for(size_t j = 0; j < n_fixed; j++)
		{	double value = sample_cov(i, j);
			double check = full_cov(i, j);
			max_err = std::max(max_err, std::fabs(value - check) / max_cov );
		}
	}
	ok &= max_err < .05;
	//
	if( ! ok )
		std::cout << "\nrandom_seed = " << random_seed << "\n";
	//
	return ok;
}
// END C++
