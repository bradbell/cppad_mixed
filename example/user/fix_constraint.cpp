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
$begin fix_constraint.cpp$$
$spell
	CppAD
	cppad
	hes
	eval
	interp
	xam
$$

$section Using Constraints: Example and Test$$

$head Model$$
$latex \[
	\B{p}( y_i | \theta , u ) \sim \B{N} ( u_i + \theta_i , 1 )
\] $$
$latex \[
	\B{p}( u_i | \theta ) \sim \B{N} ( 0 , 1 )
\] $$
$latex \[
	\B{p}( \theta ) \sim \B{U} ( - \infty , + \infty )
\] $$
where $latex \B{U} ( - \infty ,  + \infty )$$ is the improper uniform prior
on $latex [- \infty , + \infty ]$$.
It follows that the Laplace approximation is exact and
$latex \[
	\B{p}( y_i | \theta ) \sim \B{N} ( \theta_i , 2 )
\] $$
The corresponding objective for the fixed effects is equivalent to:
$latex \[
	\frac{1}{2} \sum_{i=0}^{N-1} ( y_i - \theta_i )^2
\] $$
For this problem we add the explicit constraint
$latex \[
	\frac{1}{2} \sum_i \theta_i^2 = 1;
\] $$
The corresponding Lagrangian is
$latex \[
L( \theta , \lambda ) =
	\frac{1}{2} \sum_{i=0}^{N-1} ( y_i - \theta_i )^2
		+ \lambda \left( \frac{1}{2} \sum_i \theta_i^2 - 1 \right)
\] $$

$code
$srcfile%example/user/fix_constraint.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>

namespace {
	using CppAD::log;
	using CppAD::AD;
	//
	using CppAD::mixed::sparse_mat_info;
	using CppAD::mixed::a1_double;
	using CppAD::mixed::a2_double;
	using CppAD::mixed::d_vector;
	using CppAD::mixed::a1_vector;
	using CppAD::mixed::a2_vector;
	//
	class mixed_derived : public cppad_mixed {
	private:
		size_t                n_fixed_;
		const d_vector&       y_;
	public:
		// constructor
		mixed_derived(
			size_t                   n_fixed        ,
			size_t                   n_random       ,
			bool                     quasi_fixed    ,
			bool                     bool_sparsity  ,
			const  sparse_mat_info&  A_info         ,
			const d_vector&          y              ) :
			cppad_mixed(
				n_fixed, n_random, quasi_fixed, bool_sparsity, A_info
			),
			n_fixed_(n_fixed),
			y_(y)
		{}
		// implementation of ran_likelihood
		virtual a2_vector ran_likelihood(
			const a2_vector&         theta  ,
			const a2_vector&         u      )
		{	assert( u.size() == y_.size() );
			assert( theta.size() == y_.size() );
			a2_vector vec(1);

			// initialize part of log-density that is always smooth
			vec[0] = a2_double(0.0);

			// pi
			// sqrt_2pi = CppAD::sqrt(8.0 * CppAD::atan(1.0) );

			for(size_t i = 0; i < y_.size(); i++)
			{	a2_double mu     = u[i] + theta[i];
				a2_double sigma  = a2_double(1.0);
				a2_double res    = (y_[i] - mu) / sigma;

				// p(y_i | u, theta)
				vec[0] += res*res / a2_double(2.0);
				// following term does not depend on fixed or random effects
				// vec[0] += log(sqrt_2pi * sigma);

				// p(u_i | theta)
				vec[0] += u[i] * u[i] / a2_double(2.0);
				// following term does not depend on fixed or random effects
				// vec[0] += log(sqrt_2pi);
			}
			return vec;
		}
		// ------------------------------------------------------------------
		// ran_likelihood
		// fix_constraint
		virtual a1_vector fix_constraint(
			const a1_vector&         fixed_vec  )
		{	a1_vector ret_val(1);
			//
			ret_val[0] = 0.0;
			for(size_t i = 0; i < fixed_vec.size(); i++)
				ret_val[0] += fixed_vec[i] * fixed_vec[i];
			ret_val[0] /= 2.0;
			//
			return ret_val;
		}
	};
}

bool fix_constraint_xam(void)
{
	bool   ok = true;
	double inf = std::numeric_limits<double>::infinity();
	double tol = 1e-8;

	size_t n_data   = 3;
	size_t n_fixed  = n_data;
	size_t n_random = n_data;
	d_vector
		fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
	for(size_t i = 0; i < n_fixed; i++)
	{	fixed_lower[i] = - inf;
		fixed_in[i]    = 0.1;
		fixed_upper[i] = inf;
	}
	//
	// explicit constriants (in addition to l1 terms)
	d_vector fix_constraint_lower(1), fix_constraint_upper(1);
	fix_constraint_lower[0] = 1.0;
	fix_constraint_upper[0] = 1.0;
	//
	d_vector data(n_data), random_in(n_random);
	for(size_t i = 0; i < n_data; i++)
	{	data[i]       = double(i + 1);
		random_in[i] = 0.0;
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed   = false;
	bool bool_sparsity = false;
	CppAD::mixed::sparse_mat_info A_info; // empty matrix
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_info, data
	);
	mixed_object.initialize(fixed_in, random_in);

	// optimize the fixed effects using full Newton method
	std::string fixed_ipopt_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           second-order\n"
		"String  derivative_test_print_all yes\n"
		"Numeric tol                       1e-8\n"
		"Integer max_iter                  15\n"
	;
	// random_ipopt_options is non-empty, so using ipopt for random effects
	std::string random_ipopt_options =
		"Integer print_level 0\n"
		"String  sb          yes\n"
		"String  derivative_test second-order\n"
	;
	// lower and upper limits for random effects
	d_vector random_lower(n_random), random_upper(n_random);
	for(size_t i = 0; i < n_random; i++)
	{	random_lower[i] = -inf;
		random_upper[i] = +inf;
	}
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
	d_vector fixed_out = solution.fixed_opt;
	//
	// check constraint
	double sum = 0.0;
	for(size_t i = 0; i < n_fixed; i++)
		sum += fixed_out[i] * fixed_out[i];
	ok &= fabs( sum / 2.0 - 1.0 ) <= tol;

	// compute lagranges multiplier by averaging
	sum = 0.0;
	for(size_t i = 0; i < n_fixed; i++)
		sum += (fixed_out[i] - data[i]) / fixed_out[i];
	double lambda = sum / double(n_fixed);

	// check partials of Lagragian w.r.t fixed effects
	for(size_t i = 0; i < n_fixed; i++)
	{	double err	= data[i] - fixed_out[i] + lambda * fixed_out[i];
		ok         &= fabs(err) < tol;
	}
	return ok;
}
// END C++
