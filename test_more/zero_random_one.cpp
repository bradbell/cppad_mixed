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
$begin zero_random_one.cpp$$
$spell
	CppAD
	cppad
	hes
	eval
	interp
	xam
$$

$section Case Where Optimum is Non-Zero Fixed and Zero Random Effects$$

$head Model$$
$latex \[
	\B{p}( y_0 | \theta , u ) \sim \B{N} left( \theta , \sigma_y^2 \right)
\] $$
$latex \[
	\B{p}( y_1 | \theta , u ) \sim \B{N} left( u + \theta , \sigma_y^2 \right)
\] $$
$latex \[
	\B{p}( u | \theta ) \sim \B{N} \left( 0 , \sigma_u^2 \right)
\] $$
It follows that the Laplace approximation is exact and
$latex \[
	\B{p}( y_0 | \theta ) \sim \B{N} \left( \theta , sigma_y^2 \right)
\] $$
$latex \[
\B{p}( y_1 | \theta ) \sim \B{N} \left( \theta , \sigma_u^2 + \sigma_y^2 \right)
\] $$
The corresponding objective for the fixed effects is equivalent to:
$latex \[
F( \theta ) = \frac{1}{2} \left[
	\sigma_y^{-2} ( y_0 - \theta )^2 +
	\left( \sigma_y^2 + \sigma_u^2 \right)^{-1} ( y_1 - \theta )^2
\right]
\] $$
The constraints on the fixed effect are
$latex \[
	- \infty \leq \theta \leq + \infty
\] $$

$code
$srcfile%test_more/zero_random_one.cpp
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
	using CppAD::mixed::sparse_rcv;
	//
	using CppAD::mixed::a1_double;
	using CppAD::mixed::a2_double;

	class mixed_derived : public cppad_mixed {
	private:
		const size_t n_fixed_;
		const size_t n_random_;
		const double sigma_u_;
		const double sigma_y_;
		const vector<double> y_;
	// ----------------------------------------------------------------------
	public:
		// constructor
		mixed_derived(
			size_t n_fixed                    ,
			size_t n_random                   ,
			const CppAD::mixed::sparse_rcv&      A_rcv,
			double sigma_u                    ,
			double sigma_y                    ,
			const vector<double>& y           ) :
			// quasi_fixed = false, bool_sparsity = true
			cppad_mixed(n_fixed, n_random, false, true, A_rcv) ,
			n_fixed_(n_fixed)                          ,
			n_random_(n_random)                        ,
			sigma_u_(sigma_u)                          ,
			sigma_y_(sigma_y)                          ,
			y_(y)
		{	assert( n_fixed_ == 1);
			assert( n_random_ == 1 );
			assert( y_.size() == 2 );
		}
	// ----------------------------------------------------------------------
		// implementation of ran_likelihood
		template <typename Vector>
		Vector template_ran_likelihood(
			const Vector& theta  ,
			const Vector& u      )
		{	typedef typename Vector::value_type scalar;

			assert( theta.size() == n_fixed_ );
			assert( u.size()     == n_random_ );
			Vector vec(1);

			// initialize part of log-density that is always smooth
			vec[0] = scalar(0.0);

			// square root of 2 * pi
			scalar sqrt_2pi = scalar(
				 CppAD::sqrt(8.0 * CppAD::atan(1.0)
			));

			// p(y_1 | theta, u)
			scalar mu    = u[0] + theta[0];
			scalar res   = (y_[1] - mu) / sigma_y_;
			vec[0]     += log(sqrt_2pi * sigma_y_) + res * res / scalar(2.0);

			// p(u | theta)
			res     = u[0] / sigma_u_;
			vec[0] += log(sqrt_2pi) + res * res / scalar(2.0);

			return vec;
		}
		// a2_vector version of ran_likelihood
		// a3_vector version of ran_likelihood
		virtual a3_vector ran_likelihood(
			const a3_vector& fixed_vec, const a3_vector& random_vec
		)
		{	return template_ran_likelihood( fixed_vec, random_vec ); }
		// implementation of fix_likelihood
		template <typename Vector>
		Vector template_fix_likelihood(
			const Vector& fixed_vec  )
		{	typedef typename Vector::value_type scalar;

			Vector vec(1);

			// initialize part of log-density that is smooth
			vec[0] = scalar(0.0);

			// compute these factors once
			scalar sqrt_2pi = scalar(
				 CppAD::sqrt( 8.0 * CppAD::atan(1.0)
			));

			// p(y_0 | theta)
			scalar mu     = fixed_vec[0];
			scalar res    = (y_[0] - mu) / sigma_y_;
			vec[0] += log(sqrt_2pi * sigma_y_) + res * res / scalar(2.0);

			return vec;
		}
		// a1_vector version of fix_likelihood
		virtual a1_vector fix_likelihood(const a1_vector& fixed_vec)
		{	return template_fix_likelihood( fixed_vec ); }
	// ----------------------------------------------------------------------
	public:
		// User defined virtual functions
		//
		// ------------------------------------------------------------------
	};
}

bool zero_random_one(void)
{
	bool   ok = true;
	double inf = std::numeric_limits<double>::infinity();
	double tol = 1e-8;

	double sigma_u  = 0.1;
	double sigma_y  = 0.1 * 0.05;
	size_t n_fixed  = 1;
	size_t n_random = 1;;
	size_t n_data   = 2;
	vector<double>
		fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
	fixed_lower[0] = - inf; fixed_in[0] = 0.5;  fixed_upper[0] = inf;
	//
	// explicit constriants (in addition to l1 terms)
	vector<double> fix_constraint_lower(0), fix_constraint_upper(0);
	//
	vector<double> data(n_data);
	data[0] = 0.05;
	data[1] = 0.05;
	//
	vector<double>	random_in(n_random);
	random_in[0] = 0.1;

	// object that is derived from cppad_mixed
	CppAD::mixed::sparse_rcv A_rcv; // empty matrix
	mixed_derived mixed_object(
		n_fixed, n_random, A_rcv, sigma_u, sigma_y, data
	);
	mixed_object.initialize(fixed_in, random_in);

	// optimize the fixed effects
	std::string fixed_ipopt_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           adaptive\n"
		"String  derivative_test_print_all yes\n"
		"Numeric tol                       1e-8\n"
		"Integer max_iter                  15\n"
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
	vector<double> fixed_scale = fixed_in;
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
	vector<double> fixed_out = solution.fixed_opt;
	//
	// results of optimization
	double theta = fixed_out[0];

	// compute partials of F
	double F_theta = (theta - data[0]) / (sigma_y * sigma_y);
	F_theta       += (theta - data[1]) / (sigma_y*sigma_y + sigma_u*sigma_u);

	// Note that no constraints are active, (not even the l1 terms)
	// so the partials should be zero.
	ok &= fabs( F_theta ) <= 2.0 * tol;

	// Compute the optimal random effects
	vector<double> random_out = mixed_object.optimize_random(
		random_ipopt_options, fixed_out, random_lower, random_upper, random_in
	);
	double u = random_out[0];

	// partial of p(u | theta) w.r.t u
	double p_u = (theta + u - data[0]) / (sigma_y * sigma_y);
	ok &= fabs( u )   <= 2.0 * tol;
	ok &= fabs( p_u ) <= 2.0 * tol;

	return ok;
}
// END C++
