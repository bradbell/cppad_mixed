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
$begin optimize_fixed.cpp$$
$spell
	CppAD
	cppad
	hes
	eval
	interp
	xam
$$

$section Optimize Fixed Effects: Example and Test$$

$head Model$$
$latex \[
	\B{p}( y_i | \theta , u ) \sim \B{N} ( u_i + \theta_0 , \theta_1^2 )
\] $$
$latex \[
	\B{p}( u_i | \theta ) \sim \B{N} ( 0 , 1 )
\] $$
$latex \[
	\B{p}( \theta ) \sim \B{N} ( 4 , 1 )
\] $$
It follows that the Laplace approximation is exact and
$latex \[
	\B{p}( y_i | \theta ) \sim \B{N} \left( \theta_0 , 1 + \theta_1^2 \right)
\] $$
The constraints on the fixed effect are
$latex \[
	- \infty \leq \theta_0 \leq + \infty
	\R{\; and \;}
	0.1 \leq \theta_1 \leq 100
\] $$

$head Objective$$
The corresponding objective for the fixed effects is equivalent to:
$latex \[
F( \theta ) = \frac{1}{2} \left[
	( \theta_0 - 4 )^2 + ( \theta_1 - 4 )^2 +
		N \log \left( 1 + \theta_1^2 \right) +
		( 1 + \theta_1^2)^{-1} \sum_{i=0}^{N-1} ( y_i - \theta_0 )^2
\right]
\] $$

$head First Order Partials$$
The first order partial derivatives of the objective are:
$latex \[
F_0 ( \theta )
=
( \theta_0 - 4 )
-
( 1 + \theta_1^2)^{-1} \sum_{i=0}^{N-1} ( y_i - \theta_0 )
\] $$
$latex \[
F_1 ( \theta )
=
( \theta_1 - 4 )
+
N \left( 1 + \theta_1^2 \right)^{-1} \theta_1
-
( 1 + \theta_1^2)^{-2} \theta_1 \sum_{i=0}^{N-1} ( y_i - \theta_0 )^2
\] $$

$head Source Code$$
$code
$srcfile%example/user/optimize_fixed.cpp
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
	//
	typedef AD<double>    a1_double;
	typedef AD<a1_double> a2_double;

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
		// implementation of ran_likelihood
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

			for(size_t j = 0; j < n_fixed_; j++)
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

bool optimize_fixed_xam(void)
{
	bool   ok = true;
	double inf = std::numeric_limits<double>::infinity();
	double tol = 1e-8;

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
	std::string fixed_ipopt_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           adaptive\n"
		"String  derivative_test_print_all yes\n"
		"Numeric tol                       1e-8\n"
	;
	std::string random_ipopt_options =
		"Integer print_level     0\n"
		"String  sb              yes\n"
		"String  derivative_test second-order\n"
		"Numeric tol             1e-8\n"
	;
	//
	vector<double> random_lower(n_random), random_upper(n_random);
	for(size_t i = 0; i < n_random; i++)
	{	random_lower[i] = -inf;
		random_upper[i] = +inf;
	}
	// ------------------------------------------------------------------
	// use ipopt for random effects optimization
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

	// results of optimization
	double theta_0 = fixed_out[0];
	double theta_1 = fixed_out[1];

	// compute partials of F
	double sum   = 0.0;
	double sumsq = 0.0;
	for(size_t i = 0; i < n_data; i++)
	{	double diff = data[i] - theta_0;
		sum   += diff;
		sumsq += diff * diff;
	}
	double den = 1.0 + theta_1 * theta_1;
	double F_0 = (theta_0 - 4.0) - sum / den;
	double F_1 = theta_1 - 4.0;
	F_1       += double(n_data) * theta_1 / den;
	F_1       -= theta_1 * sumsq / (den * den);

	// Note that no constraints are active, (not even the l1 terms)
	// so the partials should be zero.
	ok &= CppAD::abs( F_0 ) <= tol;
	ok &= CppAD::abs( F_1 ) <= tol;
	//
	return ok;
}
// END C++
