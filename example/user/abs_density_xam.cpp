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
$begin abs_density_xam.cpp$$
$spell
	CppAD
	cppad
	hes
	eval
	interp
	xam
$$

$section Absolute Value In Log-Density: Example and Test$$.

$head Model$$
$latex \[
	\B{p}( z_i | \theta ) \sim \B{L} ( \theta_i , \sigma )
\] $$
where $latex \B{L} ( \mu , \sigma )$$ is the Laplace distribution
with mean $latex \mu$$ and standard deviation $latex \sigma$$.
The corresponding fixed likelihood
$cref/g(theta)/theory/Fixed Likelihood, g(theta)/$$
is
$latex \[
g( \theta ) = \sum_{i} \left[
	\log ( \sigma \sqrt{2} )
	+
	\sqrt{2} \; \left| \frac{ z_i - \exp( \theta_i )}{\sigma} \right|
\right]
\] $$
The optimal solution, with no constraints and no prior on $latex \theta$$ is
$latex \[
	\hat{\theta}_i = \log( z_i )
\] $$

$code
$srcfile%example/user/abs_density_xam.cpp
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
		size_t                n_fixed_;
		double                sigma_;
		const vector<double>& z_;
	public:
		// constructor
		mixed_derived(
			size_t n_fixed                    ,
			size_t n_random                   ,
			bool   quasi_fixed                ,
			const  sparse_mat_info& A_info    ,
			double sigma                      ,
			const vector<double>& z           ) :
			cppad_mixed(n_fixed, n_random, quasi_fixed, A_info) ,
			n_fixed_(n_fixed)                           ,
			sigma_(sigma)                               ,
			z_(z)
		{	assert(z.size() == n_fixed); }
	private:
		// implementation of fix_likelihood as p(z|theta) * p(theta)
		template <class Float>
		vector<Float> implement_fix_likelihood(
			const vector<Float>& fixed_vec  )
		{
			// initialize log-density
			vector<Float> vec(1 + n_fixed_);
			vec[0] = Float(0.0);

			// compute this factors once
			Float sqrt_2 = Float( CppAD::sqrt( 2.0 ) );

			for(size_t j = 0; j < n_fixed_; j++)
			{	// Data term
				Float res   = z_[j] - CppAD::exp( fixed_vec[j] );
				res        /= Float( sigma_ );
				vec[0]     += log(sigma_ * sqrt_2);
				vec[1 + j] += sqrt_2 * res;
			}
			return vec;
		}
	public:
		// ------------------------------------------------------------------
		// User defined virtual functions
		virtual vector<a1_double> fix_likelihood(
			const vector<a1_double>& fixed_vec  )
		{	return implement_fix_likelihood(fixed_vec); }
		// ------------------------------------------------------------------
	};
}

bool abs_density_xam(void)
{
	bool   ok = true;
	double inf = std::numeric_limits<double>::infinity();
	double tol = 1e-8;

	// fixed effects
	size_t n_fixed  = 3;
	vector<double>
		fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
	for(size_t j = 0; j < n_fixed; j++)
	{	fixed_lower[j] = - inf;
		fixed_in[j]    = 0.0;
		fixed_upper[j] = inf;
	}
	//
	// no random effects
	size_t n_random = 0;
	vector<double> random_in(0);
	vector<double> random_lower(n_random), random_upper(n_random);
	std::string random_ipopt_options = "";
	CppAD::mixed::box_newton_options random_box_options;
	CppAD::mixed::sparse_mat_info A_info; // empty matrix
	//
	// no constriants
	vector<double> fix_constraint_lower(0), fix_constraint_upper(0);
	//
	vector<double> z(n_fixed);
	for(size_t i = 0; i < n_fixed; i++)
		z[i] = double(i+3);

	// object that is derived from cppad_mixed
	bool quasi_fixed = false;
	double sigma     = 1.0;
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, A_info, sigma, z
	);
	mixed_object.initialize(fixed_in, random_in);

	// optimize the fixed effects using quasi-Newton method
	std::string fixed_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           second-order\n"
		"String  derivative_test_print_all yes\n"
		"Numeric tol                       1e-8\n"
	;
	CppAD::mixed::fixed_solution solution = mixed_object.optimize_fixed(
		fixed_options,
		random_box_options,
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

	for(size_t j = 0; j < n_fixed; j++)
		ok &= CppAD::abs( fixed_out[j] - CppAD::log( z[j] ) ) <= tol;

	return ok;
}
// END C++
