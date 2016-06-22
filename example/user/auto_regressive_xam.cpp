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
$begin auto_regressive_xam.cpp$$
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
$srcfile%example/user/sample_fixed_xam.cpp
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
			{	a2_double mu     = u[i] + theta[0] * double(i + 1);
				a2_double sigma  = theta[1];
				a2_double res    = (y_[i] - mu) / sigma;

				// p(y_i | u, theta)
				vec[0] += log(sigma) + res * res / a2_double(2.0);
				// following term does not depend on fixed or random effects
				// vec[0] += log(sqrt_2pi);

				// p(u_i | theta)
				if( i == 0 )
					vec[0] += u[i] * u[i] / a2_double(2.0);
				else
				{	res     = u[i] - u[i-1];
					vec[0] += res * res / a2_double(2.0);
				}
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

bool auto_regressive_xam(void)
{
	bool   ok = true;
	double inf = std::numeric_limits<double>::infinity();
	//
	// initialize gsl random number generator
	size_t random_seed = CppAD::mixed::new_gsl_rng(0);
	//
	size_t n_data   = 100;
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
	{	data[i]      = double(i + 1);
		random_in[i] = 0.0;
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed = true;
	CppAD::mixed::sparse_mat_info A_info; // empty matrix
	mixed_derived mixed_object(n_fixed, n_random, quasi_fixed, A_info, data);
	double start_seconds = CppAD::elapsed_seconds();
	std::map<std::string, size_t> size_map =
		mixed_object.initialize(fixed_in, random_in);
	double end_seconds = CppAD::elapsed_seconds();
	// std::cout << "initilaize_time = " << end_seconds - start_seconds << "\n";
	size_t before = size_map["num_bytes_before"];
	size_t after  = size_map["num_bytes_after"];
	// std::cout << "initilaize_bytes = " << before - after << std::endl;

	// optimize the fixed effects using quasi-Newton method
	std::string fixed_ipopt_options =
		"Integer print_level               5\n"
		"String  sb                        yes\n"
		"String  derivative_test           none\n"
		"Numeric tol                       1e-7\n"
	;
	std::string random_ipopt_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           second-order\n"
		"Numeric tol                       1e-7\n"
	;
	vector<double> random_lower(n_random), random_upper(n_random);
	for(size_t i = 0; i < n_random; i++)
	{	random_lower[i] = -inf;
		random_upper[i] = +inf;
	}
	// optimize fixed effects
	start_seconds = CppAD::elapsed_seconds();
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
	end_seconds = CppAD::elapsed_seconds();
	// std::cout << "optimize_fixed_seconds = "
	//	<< end_seconds - start_seconds << "\n";
	//
	// check that none of the constraints are active
	// (Note that the Lagragian w.r.t. theta[2] will be zero because
	// it does not affect the objective).
	ok &= solution.fixed_lag.size() == n_fixed;
	for(size_t i = 0; i < n_fixed; i++)
		ok &= solution.fixed_lag[i] == 0.0;
	ok &= solution.fix_con_lag.size() == 0;
	ok &= solution.ran_con_lag.size() == 0.0;
	ok &= CppAD::abs( solution.fixed_opt[0] - 1.0 ) < 1e-1;
	ok &= CppAD::abs( solution.fixed_opt[1] ) < 1e-1;
	//
	std::cout << "solution.fixed_opt = " << solution.fixed_opt << "\n";
	//
	if( ! ok )
		std::cout << "\nrandom_seed = " << random_seed << "\n";
	//
	CppAD::mixed::free_gsl_rng();
	return ok;
}
// END C++
