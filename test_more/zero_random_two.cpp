// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin zero_random_two.cpp$$
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
	\B{p}( y_i | \theta , u ) \sim \B{N} left( \theta_0 , \theta_1^2 \right)
\] $$
$latex \[
	\B{p}( u_i | \theta ) \sim \B{N} ( 0 , 1 )
\] $$
$latex \[
	\B{p}( \theta ) \sim \B{N} ( 1 , 1 )
\] $$
It follows that the Laplace approximation is exact and
$latex \[
	\B{p}( y_i | \theta ) \sim \B{N} \left( \theta_0 , \theta_1^2 \right)
\] $$
The corresponding objective for the fixed effects is equivalent to:
$latex \[
F( \theta ) = \frac{1}{2} \left[
	( \theta_0 - 1 )^2 + ( \theta_1 - 1 )^2 +
		N \log \left( \theta_1^2 \right) +
		\theta_1^{-2} \sum_{i=0}^{N-1} ( y_i - \theta_0 )^2
\right]
\] $$
The constraints on the fixed effect are
$latex \[
	- \infty \leq \theta_0 \leq + \infty
	\R{\; and \;}
	0.1 \leq \theta_1 \leq 100
\] $$

$code
$verbatim%test_more/zero_random_two.cpp
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

	class mixed_derived : public cppad_mixed {
	private:
		const size_t n_fixed_;
		const size_t n_random_;
		const vector<double> y_;
	// ----------------------------------------------------------------------
	public:
		// constructor
		mixed_derived(
			size_t n_fixed                    ,
			size_t n_random                   ,
			const vector<double>& y           ) :
			// quasi_fixed = false
			cppad_mixed(n_fixed, n_random, false) ,
			n_fixed_(n_fixed)                          ,
			n_random_(n_random)                        ,
			y_(y)
		{	assert( n_fixed_ == 2);}
	// ----------------------------------------------------------------------
	private:
		// implementation of ran_likelihood
		template <class Float>
		vector<Float> implement_ran_likelihood(
			const vector<Float>& theta  ,
			const vector<Float>& u      )
		{	assert( theta.size() == n_fixed_ );
			assert( u.size()     == n_random_ );
			vector<Float> vec(1);

			// initialize part of log-density that is always smooth
			vec[0] = Float(0.0);

			// square root of 2 * pi
			Float sqrt_2pi = Float( CppAD::sqrt(8.0 * CppAD::atan(1.0) ) );

			for(size_t i = 0; i < u.size(); i++)
			{	// p(u_i | theta)
				vec[0] += log(sqrt_2pi) + u[i] * u[i] / Float(2.0);
			}
			return vec;
		}
		// implementation of fix_likelihood
		template <class Float>
		vector<Float> implement_fix_likelihood(
			const vector<Float>& fixed_vec  )
		{	vector<Float> vec(1);

			// initialize part of log-density that is smooth
			vec[0] = Float(0.0);

			// compute these factors once
			Float sqrt_2pi = Float( CppAD::sqrt( 8.0 * CppAD::atan(1.0) ) );

			// prior for fixed effects
			for(size_t j = 0; j < n_fixed_; j++)
			{	Float mu     = Float(1.0);
				Float sigma  = Float(1.0);
				Float res    = (fixed_vec[j] - mu) / sigma;

				// p(theta_j )
				vec[0]  += log(sqrt_2pi * sigma) + res * res / Float(2.0);
			}


			// data given fixed effects
			for(size_t i = 0; i < y_.size(); i++)
			{	Float mu     = fixed_vec[0];
				Float sigma  = fixed_vec[1];
				Float res    = (y_[i] - mu) / sigma;

				// p(y_i | theta)
				vec[0] += log(sqrt_2pi * sigma) + res * res / Float(2.0);
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
		//
		virtual vector<a1_double> fix_constraint(
			const vector<a1_double>& fixed_vec  )
		{	return vector<a1_double>(0); } // empty vector
		//
		virtual void fatal_error(const std::string& error_message)
		{	std::cerr << "Error: " << error_message << std::endl;
			std::exit(1);
		}
		//
		virtual void warning(const std::string& warning_message)
		{	std::cerr << "Warning: " << warning_message << std::endl;
		}
		// ------------------------------------------------------------------
	};
}

bool zero_random_two(void)
{
	bool   ok = true;
	double inf = std::numeric_limits<double>::infinity();
	double tol = 1e-8;

	size_t n_data   = 10;
	size_t n_fixed  = 2;
	size_t n_random = 3;;
	vector<double>
		fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
	fixed_lower[0] = - inf; fixed_in[0] = 2.0; fixed_upper[0] = inf;
	fixed_lower[1] = .01;   fixed_in[1] = 0.5; fixed_upper[1] = inf;
	//
	// explicit constriants (in addition to l1 terms)
	vector<double> fix_constraint_lower(0), fix_constraint_upper(0);
	//
	vector<double> data(n_data);
	for(size_t i = 0; i < n_data; i++)
		data[i]       = double(i + 1);
	//
	vector<double>	random_in(n_random);
	for(size_t i = 0; i < n_random; i++)
		random_in[i] = 1.0;

	// object that is derived from cppad_mixed
	mixed_derived mixed_object(n_fixed, n_random, data);
	mixed_object.initialize(fixed_in, random_in);

	// optimize the fixed effects
	std::string fixed_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           second-order\n"
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
	vector<double> fixed_out = mixed_object.optimize_fixed(
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

	// results of optimization
	double theta_0 = fixed_out[0];
	double theta_1 = fixed_out[1];

	// compute partials of F
	double sum   = 0.0;
	double sumsq = 0.0;
	for(size_t i = 0; i < n_data; i++)
	{	sum   += theta_0 - data[i];
		sumsq += (theta_0 - data[i]) * (theta_0 - data[i]);
	}
	double F_0 = (theta_0 - 1.0) + sum / (theta_1 * theta_1);
	double F_1 = theta_1 - 1.0;
	F_1       += double(n_data) / theta_1;
	F_1       -= sumsq  / (theta_1 * theta_1 * theta_1);

	// Note that no constraints are active, (not even the l1 terms)
	// so the partials should be zero.
	ok &= CppAD::abs( F_0 ) <= 2.0 * tol;
	ok &= CppAD::abs( F_1 ) <= 2.0 * tol;

	// Compute the optimal random effects
	vector<double> random_out = mixed_object.optimize_random(
		random_options, fixed_out, random_lower, random_upper, random_in
	);
	// partial of p(u | theta) w.r.t u_i is equal to u_i
	for(size_t i = 0; i < n_random; i++)
		ok &= CppAD::abs( random_out[i] ) <= 2.0 * tol;

	return ok;
}
// END C++
