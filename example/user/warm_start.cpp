/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin warm_start.cpp$$
$spell
	Optimizer
	vec
$$

$nospell
$bold This is cppad_mixed--20220519 documentation:$$ Here is a link to its
$href%https://cppad-mixed.readthedocs.io%current documentation%$$.
$$
$section Warm Starting Optimization: Example and Test$$.

$head Model$$
$latex \[
	\B{p}( z_i | \theta ) \sim \B{N} ( \theta_i , 1 )
\] $$
with no prior on $latex \theta$$.
The corresponding fixed likelihood
$cref/g(theta)/theory/Fixed Likelihood, g(theta)/$$
is
$latex \[
g( \theta ) = \frac{1}{2} \sum_{i} \left[
	\log ( 2 \pi ) + ( z_i - \theta_i )^2
\right]
\] $$
We do not include the constant term $latex \log( 2 \pi )$$
in the fixed likelihood.
The optimal solution (with no constraints) is
$latex \[
	\hat{\theta}_i = z_i
\] $$

$head Bounds$$
We add lower and upper bounds that are not active at the optimal solution.
To be specific
$latex \[
	0 \leq \theta_i \leq z_i + 1
\] $$

$head Maximum Iterations$$
We use 5 for the maximum number of iterations so that the optimization
problem does not solve on the first try.
A warm start is used and the problem does solve within the limit
of another 5 iterations.

$head Optimizer Trace$$
This example uses the optimizer trace information; see
$cref/trace_vec/fixed_solution/trace_vec/$$.

$code
$srcthisfile%0%// BEGIN C++%// END C++%1%$$
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
	using CppAD::mixed::d_sparse_rcv;
	using CppAD::mixed::a1_double;
	using CppAD::mixed::d_vector;
	using CppAD::mixed::a1_vector;
	//
	class mixed_derived : public cppad_mixed {
	private:
		const size_t          n_fixed_;
		const d_vector&       z_;
	public:
		// constructor
		mixed_derived(
			size_t                 n_fixed        ,
			size_t                 n_random       ,
			const d_vector&        z              ) :
			cppad_mixed(n_fixed, n_random)  ,
			n_fixed_(n_fixed)               ,
			z_(z)
		{	assert(z.size() == n_fixed); }
		//
		// implementation of fix_likelihood as p(z|theta)
		a1_vector fix_likelihood(
			const a1_vector&         fixed_vec  ) override
		{
			// initialize log-density
			a1_vector vec(1);
			vec[0] = 0.0;

			for(size_t j = 0; j < n_fixed_; j++)
			{
				// Data term p(z|theta)
				a1_double res  = (z_[j] - fixed_vec[j]);
				vec[0]    += res * res / 2.0;
			}
			return vec;
		}
		//
		// warning
		bool   suppress_warning_;
		size_t warning_count_;
		void warning(const std::string& warning_message) override
		{	++warning_count_;
			if( ! suppress_warning_ )
			std::cerr << "cppad_mixed warning: " << warning_message << "\n";
		}
	};
}

bool warm_start_xam(void)
{
	bool   ok = true;
	double inf = std::numeric_limits<double>::infinity();
	double tol = 1e-8;
	//
	// n_fixed
	size_t n_fixed  = 3;
	//
	// z
	d_vector z(n_fixed);
	for(size_t j = 0; j < n_fixed; j++)
		z[j] = double(j+1);
	//
	// fixed_lower, fixed_in, fixed_upper
	d_vector fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
	for(size_t j = 0; j < n_fixed; j++)
	{	fixed_lower[j] = 0.0;
		fixed_in[j]    = 0.0;
		fixed_upper[j] = z[j] + 1.0;
	}
	//
	// n_random, random_in
	size_t n_random = 0;
	d_vector random_in(0);
	//
	// fix_constraint_lower, fix_constraint_upper
	d_vector fix_constraint_lower(0), fix_constraint_upper(0);
	//
	// mixed_object
	mixed_derived mixed_object(n_fixed, n_random, z);
	mixed_object.initialize(fixed_in, random_in);
	//
	// optimize the fixed effects using quasi-Newton method
	std::string fixed_ipopt_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           first-order\n"
		"String  derivative_test_print_all yes\n"
		"Numeric tol                       1e-8\n"
		"Integer max_iter                  5\n"
	;
	std::string random_ipopt_options =
		"Integer print_level 0\n"
		"String  sb          yes\n"
		"String  derivative_test second-order\n"
	;
	//
	// random_lower, random_upper
	d_vector random_lower(n_random), random_upper(n_random);
	for(size_t i = 0; i < n_random; i++)
	{	random_lower[i] = -inf;
		random_upper[i] = +inf;
	}
	// fixed_in
	d_vector fixed_scale = fixed_in;
	//
	// first optimization attempt (max_iter not large enough)
	mixed_object.warning_count_    = 0;
	mixed_object.suppress_warning_ = true;
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
	// optimization did not converge which causes warnings
	ok &= mixed_object.warning_count_ > 0;
	//
	// should have reached the maximum number of iterations
	// (the trace includes iteraiton zero as the starting point)
	ok &= solution.trace_vec.size() == 6;
	ok &= solution.trace_vec[5].iter == 5;
	//
	// second optimzation attempt (max_iter large enough with warm start)
	mixed_object.warning_count_    = 0;
	mixed_object.suppress_warning_ = false;
	solution = mixed_object.optimize_fixed(
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
		random_in,
		solution.warm_start
	);
	// should not be any warnings this time
	ok &= mixed_object.warning_count_ == 0;
	//
	// this time optimization should have completed in 3 iterations
	ok &= solution.trace_vec.size() == 4;
	ok &= solution.trace_vec[3].iter == 3;
	//
	// final solution
	d_vector fixed_out = solution.fixed_opt;
	//
	for(size_t j = 0; j < n_fixed; j++)
		ok &= fabs( fixed_out[j] - z[j] ) <= tol;
	//
	return ok;
}
// END C++
