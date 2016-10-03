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
$begin lasso.cpp$$
$spell
	CppAD
	cppad
	hes
	eval
	interp
	xam
$$

$section Lasso on Fixed Effects: Example and Test$$.

$head Model$$
We are given a set of times
$latex \{ t_i \W{:} i = 0 , \ldots , N-1 \}$$ and
$latex \[
\begin{array}{rcl}
	q( \theta, s )
	& = &
	\theta_0 s + \theta_1 \sin ( 2 \pi s  ) + \theta_2 \cos ( 2 \pi s )
	\\
	z_i & = & q( \theta , t_i ) + e_i
	\\
	\B{p} ( e_i | \theta ) & \sim & \B{N} ( 0, \sigma )
\end{array}
\] $$
The idea in Lasso is that one or more of the components of
$latex \theta$$ are zero and using the Laplace prior
we can recover this fact.
We use $latex \B{L} ( \mu , \sigma )$$ to denote the Laplace distribution
with mean $latex \mu$$ and standard deviation $latex \sigma$$.
$latex \[
	\B{p} ( \theta ) \sim \B{L} ( 0 , \delta )
\] $$
The corresponding fixed likelihood
$cref/g(theta)/theory/Fixed Likelihood, g(theta)/$$
is
$latex \[
g( \theta ) =
\sum_{i=0}^{N-1} \left[
	\log ( \sigma \sqrt{2 \pi} )
	+
	\left( \frac{ z_i - q( \theta , t_i ) }{2 \sigma} \right)^2
\right]
+
\sum_{j=0}^2 \left[
	\log \left( \delta \sqrt{2} \right)
	+
	\sqrt{2} \; \left| \frac{\theta_j}{\delta} \right|
\right]
\] $$

$code
$srcfile%example/user/lasso.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/manage_gsl_rng.hpp>
# include <gsl/gsl_randist.h>

namespace {
	using CppAD::vector;
	using CppAD::log;
	using CppAD::AD;
	using CppAD::mixed::sparse_mat_info;
	//
	typedef AD<double>    a1_double;
	typedef AD<a1_double> a2_double;

	class mixed_derived : public cppad_mixed {
	private:
		size_t                n_fixed_;
		double                sigma_;
		double                delta_;
		const vector<double>& t_;
		const vector<double>& z_;
	public:
		// constructor
		mixed_derived(
			size_t n_fixed                    ,
			size_t n_random                   ,
			bool   quasi_fixed                ,
			const  sparse_mat_info& A_info    ,
			double sigma                      ,
			double delta                      ,
			const vector<double>& t           ,
			const vector<double>& z           ) :
			cppad_mixed(n_fixed, n_random, quasi_fixed, A_info) ,
			n_fixed_(n_fixed)                           ,
			sigma_(sigma)                               ,
			delta_(delta)                               ,
			t_(t)                                       ,
			z_(z)
		{	assert(n_fixed == 3);
			assert( t.size() == z.size() );
		}
		// implementation of fix_likelihood as p(z|theta) * p(theta)
		virtual vector<a1_double> fix_likelihood(
			const vector<a1_double>& fixed_vec  )
		{	size_t N = t_.size();

			// initialize log-density
			vector<a1_double> vec(1 + n_fixed_);
			vec[0] = a1_double(0.0);

			// compute this factors once
			a1_double   pi     = a1_double( 4.0 * CppAD::atan(1.0) );
			a1_double sqrt_2   = a1_double( CppAD::sqrt( 2.0 ) );
			// a1_double sqrt_2pi = CppAD::sqrt( 2.0 * pi );

			// Data terms p(z|theta)
			for(size_t i = 0; i < N; i++)
			{	a1_double q_i   = fixed_vec[0] * t_[i];
			    q_i        += fixed_vec[1] * sin( 2.0 * pi * t_[i] );
			    q_i        += fixed_vec[2] * cos( 2.0 * pi * t_[i] );
				a1_double res   = z_[i] - q_i;
				res        /= a1_double( sigma_ );
				vec[0]     += res * res / 2.0;
				// following term does not depend on the fixed effects
				// vec[0]     += log(sigma_ * sqrt_2pi);
			}

			// Prior terms p(theta)
			for(size_t j = 0; j < n_fixed_; j++)
			{	// following term does not depend on the fixed effects
				// vec[0]    += log( delta_ * sqrt_2 );
				vec[1 + j] = sqrt_2 * fixed_vec[j] / delta_;
			}
			return vec;
		}
	};
}

bool lasso_xam(void)
{
	bool   ok         = true;
	double inf         = std::numeric_limits<double>::infinity();
	// size_t random_seed = CppAD::mixed::new_gsl_rng(0);
	CppAD::mixed::new_gsl_rng(0);
    gsl_rng* rng       = CppAD::mixed::get_gsl_rng();

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
	CppAD::mixed::sparse_mat_info A_info; // empty matrix
	//
	// no constriants
	vector<double> fix_constraint_lower(0), fix_constraint_upper(0);
	//
	size_t n_data = 100;
	double sigma  = 0.1;
	double pi     = 4.0 * std::atan(1.0);
	vector<double> z(n_data), t(n_data);
	for(size_t i = 0; i < n_data; i++)
	{	t[i] = double(i) / double(n_data - 1) - 0.5;
		//
		// simulation theta_0 = 0, theta_1 = 1, theta_2 = 0
		double q_i = 0.0 * t[i];
		q_i       += 1.0 * sin(2.0 * pi *t[i]);
		q_i       += 0.0 * cos(2.0 * pi *t[i]);
		double e_i = gsl_ran_gaussian(rng, sigma);
		z[i]       = q_i + e_i;
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed = false;
	double delta     = 0.002;
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, A_info, sigma, delta, t, z
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
	// coefficients that should be zero
	ok &= CppAD::abs( fixed_out[0] ) <= 1e-9;
	ok &= CppAD::abs( fixed_out[2] ) <= 1e-9;
	// non-zero coefficient has shrunk (due to prior)
	ok &= fixed_out[1] < 1.0;
	ok &= 0.75 < fixed_out[1];

	//
	CppAD::mixed::free_gsl_rng();
	return ok;
}
// END C++
