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
$begin ar1_xam.cpp$$
$escape $$
$spell
	ar1_xam
	gsl_rng
	ipopt
$$

$section A First Order Auto-Regressive Example and Speed Test$$

$head Syntax$$
$codei%build/speed/ar1_xam \
	%trace_optimize_fixed% \
	%random_seed% \
	%number_random%
%$$

$head Input$$

$subhead trace_optimize_fixed$$
This is either $code true$$ or $code false$$.
If it is true, a $icode%print_level% = 5%$$
$cref/trace/ipopt_trace/$$ of the fixed effects optimization
is included in the program output.
Otherwise the ipopt $icode print_level$$ is zero and
no such trace is printed.

$subhead random_seed$$
This is a non-negative integer equal to the
seed for the random number generator,
to be specific,
$cref/s_in/manage_gsl_rng/new_gsl_rng/s_in/$$ used during the call to
$code new_gsl_rng$$.

$subhead number_random$$
THis is a positive integer specifying the number of random effects.
This is also the number of time points (and data values) in the model.

$head Output$$

$subhead actual_seed$$
If $icode random_seed$$ is zero,
the system clock, instead of $icode random_seed$$,
is used to seed the random number generator.
The actual random seed $icode actual_seed$$ is printed
so that you can reproduce results when $icode random_seed$$ is zero.

$subhead optimize_fixed_seconds$$
Is the number of seconds used by the call to
$cref optimize_fixed$$ that is used to compute the
optimal fixed effects.

$subhead ar1_xam_ok$$
If this program passes it's correctness test,
$icode ar1_xam_ok$$ is true and the program return code is $code 0$$.
Otherwise $icode ar1_xam_ok$$ it is false and the return code is $code 1$$.

$code
$srcfile%speed/ar1_xam.cpp
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
	using std::cout;
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
		// -------------------------------------------------------------------
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
	template <class Value>
	void label_print(const char* label, const Value& value)
	{	cout << std::setw(35) << std::left << label;
		cout << " = " << value << std::endl;
	}
	void label_print(const char* label, const double& value)
	{	cout << std::setw(35) << std::left << label;
		size_t n_digits = 3 + size_t( std::log10(value) + 1e-9 );
		cout << " = " << std::setprecision(n_digits) << value << std::endl;
	}
}

int main(int argc, const char* argv[])
{
	bool   ok = true;
	double inf = std::numeric_limits<double>::infinity();
	//
	const char* arg_name[] = {
		"trace_optimize_fixed",
		"random_seed",
		"number_random"
	};
	size_t n_arg = sizeof(arg_name) / sizeof(arg_name[0]);
	//
	if( size_t(argc) != 1 + n_arg )
	{	// print usage error message
		std::cerr << "usage: " << argv[0];
		for(size_t i = 0; i < n_arg; i++)
			std::cerr << " \\ \n\t" << arg_name[i];
		std::cerr << "\n";
		std::exit(1);
	}
	//
	// get command line arguments
	assert( n_arg == 3 );
	bool   trace_optimize_fixed   = std::string( argv[1] ) == "true";
	size_t random_seed            = std::atoi( argv[2] );
	size_t number_random          = std::atoi( argv[3] );
	//
	// print the command line arugments with labels for each value
	for(size_t i = 0; i < n_arg; i++)
		label_print(arg_name[i], argv[1+i]);
	//
	// initialize gsl random number generator
	size_t actual_seed = CppAD::mixed::new_gsl_rng(random_seed);
	label_print("actual_seed", actual_seed);
	//
	size_t n_data   = number_random;
	size_t n_fixed  = 2;
	size_t n_random = number_random;
	vector<double>
		fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
	fixed_lower[0] = - inf; fixed_in[0] = 2.0; fixed_upper[0] = inf;
	fixed_lower[1] = .1;    fixed_in[1] = 0.5; fixed_upper[1] = inf;
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
	mixed_object.initialize(fixed_in, random_in);

	// optimize the fixed effects using quasi-Newton method
	std::string random_ipopt_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           second-order\n"
		"Numeric tol                       1e-7\n"
	;
	std::string fixed_ipopt_options =
		"String  sb                        yes\n"
		"String  derivative_test           none\n"
		"Numeric tol                       1e-7\n"
	;
	if( trace_optimize_fixed )
		fixed_ipopt_options += "Integer print_level               5\n";
	else
		fixed_ipopt_options += "Integer print_level               0\n";
	//
	vector<double> random_lower(n_random), random_upper(n_random);
	for(size_t i = 0; i < n_random; i++)
	{	random_lower[i] = -inf;
		random_upper[i] = +inf;
	}
	// optimize fixed effects
	double start_seconds = CppAD::elapsed_seconds();
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
	double end_seconds = CppAD::elapsed_seconds();
	if( trace_optimize_fixed )
		std::cout << std::endl;
	label_print("optimize_fixed_seconds", end_seconds - start_seconds);
	//
	// Check that only the lower limit on theta[1] is active.
	ok &= solution.fixed_lag.size() == n_fixed;
	ok &= solution.fixed_lag[0] == 0.0;
	ok &= solution.fixed_lag[1] != 0.0;
	ok &= solution.fix_con_lag.size() == 0;
	ok &= solution.ran_con_lag.size() == 0;
	ok &= CppAD::abs( solution.fixed_opt[0] - 1.0 ) < 1e-1;
	ok &= CppAD::abs( solution.fixed_opt[1] - fixed_lower[1] ) < 1e-8;
	//
	if( ok )
		label_print("ar1_xam_ok", "true");
	else
		label_print("ar1_xam_ok", "false");
	//
	CppAD::mixed::free_gsl_rng();
	if( ok )
		return 0;
	return 1;
}
// END C++
