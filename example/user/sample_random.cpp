/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin sample_random.cpp$$
$spell
	CppAD
	cppad
	hes
	eval
	interp
	xam
$$

$nospell
$bold This is cppad_mixed--20220519 documentation:$$ Here is a link to its
$href%https://cppad-mixed.readthedocs.io%current documentation%$$.
$$
$section Sample From Fixed Effects Posterior: Example and Test$$

$code
$srcthisfile%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/manage_gsl_rng.hpp>
# include <Eigen/Dense>

namespace {
	using CppAD::vector;
	using CppAD::log;
	using CppAD::AD;
	//
	using CppAD::mixed::d_sparse_rcv;
	using CppAD::mixed::d_vector;

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
			const d_vector&       y              ) :
			cppad_mixed(
				n_fixed, n_random, quasi_fixed, bool_sparsity
			)                     ,
			n_fixed_(n_fixed)     ,
			n_random_(n_random)   ,
			y_(y)
		{	assert( n_fixed == 2);
			assert( y_.size() == n_random_ );
		}
	// ----------------------------------------------------------------------
		// implementation of ran_likelihood
		a1_vector ran_likelihood(
			const a1_vector&         theta  ,
			const a1_vector&         u      ) override
		{
			assert( theta.size() == n_fixed_ );
			assert( u.size() == y_.size() );
			a1_vector vec(1);

			// initialize part of log-density that is always smooth
			vec[0] = 0.0;

			// sqrt_2pi = CppAD::sqrt(8.0 * CppAD::atan(1.0) );

			for(size_t i = 0; i < n_random_; i++)
			{	a1_double mu     = u[i] + theta[0];
				a1_double sigma  = theta[1];
				a1_double res    = (y_[i] - mu) / sigma;

				// p(y_i | u, theta)
				vec[0] += log(sigma) + res * res / 2.0;
				// following term does not depend on fixed or random effects
				// vec[0] += log(sqrt_2pi);

				// p(u_i | theta)
				a1_double sq = u[i] * u[i];
				if( i == 0 )
					sq = (u[0] + u[1]) * (u[0] + u[1]);
				vec[0] += sq / 2.0;
				// following term does not depend on fixed or random effects
				// vec[0] += log(sqrt_2pi);
			}
			return vec;
		}
	};
}

bool sample_random_xam(void)
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
	//
	d_vector data(n_data), random_in(n_random);
	for(size_t i = 0; i < n_data; i++)
	{	data[i]       = double(i + 1);
		random_in[i]    = 0.0;
	}
	d_vector fixed_vec(n_fixed);
	for(size_t i = 0; i < n_fixed; i++)
		fixed_vec[i] = double(i + 1);

	// object that is derived from cppad_mixed
	bool quasi_fixed   = true;
	bool bool_sparsity = true;
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, data
	);
	mixed_object.initialize(fixed_vec, random_in);

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
	//
	// compute the optimal random effects
	d_vector random_opt = mixed_object.optimize_random(
		random_ipopt_options, fixed_vec, random_lower, random_upper, random_in
	);
	//
	// sample from the posterior for random effects given fixed effects
	// and compute the  sample covariance matrix
	size_t n_sample = 10000;
	d_vector sample(n_sample * n_random);
	std::string error_msg = mixed_object.sample_random(
		sample,
		random_ipopt_options,
		fixed_vec,
		random_lower,
		random_upper,
		random_in
	);
	ok &= error_msg == "";
	d_vector sample_cov(n_random * n_random);
	for(size_t i = 0; i < n_random; i++)
		for(size_t j = 0; j < n_random; j++)
			sample_cov[i * n_random + j] = 0;
	for(size_t i_sample = 0; i_sample < n_sample; i_sample++)
	{	d_vector diff(n_random);
		for(size_t j = 0; j < n_random; j++)
			diff[j] = sample[i_sample * n_random + j] - random_opt[j];
		for(size_t i = 0; i < n_random; i++)
			for(size_t j = 0; j < n_random; j++)
				sample_cov[i * n_random + j] +=
					diff[i] * diff[j] / double(n_sample);
	}
	//
	// For i > 1, the terms involving u[i] is
	//	  0.5 * ( (y[i] - theta[0] - u[i])^2 / theta[1]^2 + u[i]^2  )
	// The terms involving u[0] and u[1] are
	//	  0.5 * ( (y[1] - theta[0] - u[1])^2 / theta[1]^2 + u[1]^2  )
	//	+ 0.5 * ( (y[0] - theta[0] - u[0])^2 / theta[1]^2 + (u[0]+u[1])^2  )
	double diag = 1.0 / (fixed_vec[1] * fixed_vec[1]) + 1.0;
	typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > double_mat;
	double_mat info_mat = double_mat::Zero(n_random, n_random);
	for(size_t i = 0; i < n_random; i++)
		info_mat(i, i) = diag;
	info_mat(0, 1)  = 1.0;
	info_mat(1, 0)  = 1.0;
	info_mat(1, 1) += 1.0;
	//
	double_mat cov_mat = info_mat.inverse();
	double max_err = 0.0;
	for(size_t i = 0; i < n_random; i++)
	{	for(size_t j = 0; j < n_random; j++)
		{	double value = sample_cov[i * n_random + j];
			double check = cov_mat(i, j);
			max_err = std::max(max_err, fabs(value - check) / diag);
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
