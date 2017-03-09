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
$begin sample_fixed.cpp$$
$spell
	CppAD
	cppad
	hes
	eval
	interp
	xam
$$

$section Sample From Fixed Effects Posterior: Example and Test$$

$code
$srcfile%example/user/sample_fixed.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$

$end
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
	//
	using CppAD::mixed::a1_double;
	using CppAD::mixed::a2_double;

	class mixed_derived : public cppad_mixed {
	private:
		const size_t          n_fixed_;
		const size_t          n_random_;
		const vector<double>& y_;
	// ----------------------------------------------------------------------
	public:
		// constructor
		mixed_derived(
			size_t                 n_fixed       ,
			size_t                 n_random      ,
			bool                   quasi_fixed   ,
			bool                   bool_sparsity ,
			const sparse_mat_info& A_info        ,
			const vector<double>& y              ) :
			cppad_mixed(
				n_fixed, n_random, quasi_fixed, bool_sparsity, A_info
			)                     ,
			n_fixed_(n_fixed)     ,
			n_random_(n_random)   ,
			y_(y)
		{	assert( n_fixed == 3);
			assert( y_.size() == n_random_ );
		}
	// ----------------------------------------------------------------------
		// implementation of ran_likelihood
		// Note that theta[2] is not used
		virtual vector<a2_double> ran_likelihood(
			const vector<a2_double>& theta  ,
			const vector<a2_double>& u      )
		{	assert( theta.size() == n_fixed_ );
			assert( u.size() == y_.size() );
			vector<a2_double> vec(1);

			// initialize part of log-density that is always smooth
			vec[0] = a2_double(0.0);

			// sqrt_2pi = CppAD::sqrt(8.0 * CppAD::atan(1.0) );

			for(size_t i = 0; i < n_random_; i++)
			{	a2_double mu     = u[i] + theta[0];
				a2_double sigma  = theta[1];
				a2_double res    = (y_[i] - mu) / sigma;

				// p(y_i | u, theta)
				vec[0] += log(sigma) + res * res / a2_double(2.0);
				// following term does not depend on fixed or random effects
				// vec[0] += log(sqrt_2pi);

				// p(u_i | theta)
				vec[0] += u[i] * u[i] / a2_double(2.0);
				// following term does not depend on fixed or random effects
				// vec[0] += log(sqrt_2pi);
			}
			return vec;
		}
		// implementation of fix_likelihood
		virtual vector<a1_double> fix_likelihood(
			const vector<a1_double>& fixed_vec  )
		{	assert( fixed_vec.size() == n_fixed_ );
			vector<a1_double> vec(1);

			// initialize part of log-density that is smooth
			vec[0] = a1_double(0.0);

			// sqrt_2pi = CppAD::sqrt( 8.0 * CppAD::atan(1.0) );

			// Note that theta[2] is not included
			for(size_t j = 0; j < 2; j++)
			{	a1_double mu     = a1_double(4.0);
				a1_double sigma  = a1_double(1.0);
				a1_double res    = (fixed_vec[j] - mu) / sigma;

				// This is a Gaussian term, so entire density is smooth
				vec[0]  += res * res / a1_double(2.0);
				// following term does not depend on fixed effects
				// vec[0]  += log(sqrt_2pi * sigma);
			}
			return vec;
		}
	};
}

bool sample_fixed_xam(void)
{
	bool   ok = true;
	double inf = std::numeric_limits<double>::infinity();
	double eps = 10. * std::numeric_limits<double>::epsilon();
	//
	// initialize gsl random number generator
	size_t random_seed = CppAD::mixed::new_gsl_rng(0);
	//
	size_t n_data   = 10;
	size_t n_fixed  = 3;
	size_t n_random = n_data;
	vector<double>
		fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
	fixed_lower[0] = - inf; fixed_in[0] = 2.0; fixed_upper[0] = inf;
	fixed_lower[1] = .01;   fixed_in[1] = 0.5; fixed_upper[1] = inf;
	// Set upper and lower limit for theta[2] equal so that it is bounded
	// (otherwise the implicit information matrix would be singular).
	fixed_lower[2] = 1.0;   fixed_in[2] = 1.0; fixed_upper[2] = 1.0;
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
	bool quasi_fixed   = true;
	bool bool_sparsity = true;
	CppAD::mixed::sparse_mat_info A_info; // empty matrix
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_info, data
	);
	mixed_object.initialize(fixed_in, random_in);

	// optimize the fixed effects using quasi-Newton method
	std::string fixed_ipopt_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           adaptive\n"
		"String  derivative_test_print_all yes\n"
		"Numeric tol                       1e-8\n"
		"Integer max_iter                  15\n"
	;
	std::string random_ipopt_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           second-order\n"
		"Numeric tol                       1e-8\n"
	;
	vector<double> random_lower(n_random), random_upper(n_random);
	for(size_t i = 0; i < n_random; i++)
	{	random_lower[i] = -inf;
		random_upper[i] = +inf;
	}
	// optimize fixed effects
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
	//
	// check that none of the constraints are active
	// (Note that the Lagragian w.r.t. theta[2] will be zero because
	// it does not affect the objective).
	ok &= solution.fixed_lag.size() == n_fixed;
	for(size_t i = 0; i < n_fixed; i++)
		ok &= solution.fixed_lag[i] == 0.0;
	ok &= solution.fix_con_lag.size() == 0;
	ok &= solution.ran_con_lag.size() == 0.0;
	//
	// corresponding optimal random effects
	vector<double> random_opt = mixed_object.optimize_random(
		random_ipopt_options,
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
	size_t n_sample = 20000;
	CppAD::vector<double> sample( n_sample * n_fixed );
	mixed_object.sample_fixed(
		sample,
		information_info,
		solution,
		fixed_lower,
		fixed_upper,
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
	matrix info_mat(n_fixed-1, n_fixed-1);
	// note theta[2] does not have any non-zero terms in Hessian
	size_t K = ( (n_fixed-1) * n_fixed ) / 2;
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
			if( i == 2 || j == 2 )
				ok &= std::fabs(value) < eps;
			else
			{	double check = cov_mat(i, j);
				double scale = std::sqrt( cov_mat(i, i) * cov_mat(j, j) );
				ok &= std::fabs(value - check) / scale < .05;
			}
		}
	}
	//
	if( ! ok )
		std::cout << "\nrandom_seed = " << random_seed << "\n";
	//
	CppAD::mixed::free_gsl_rng();
	return ok;
}
// END C++
