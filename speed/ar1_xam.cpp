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
	cppad
	CppAD
$$

$section A First Order Auto-Regressive Example and Speed Test$$

$head Syntax$$
$codei%build/speed/ar1_xam \
	%trace_optimize_fixed% \
	%ipopt_solve% \
	%random_seed% \
	%number_times%
%$$

$head Problem$$

$subhead Data$$
For $latex t = 0 , \ldots , T - 1$$,
$latex y_t = (1 + t) + e_t $$,
where $latex e_t \sim \B{N}( 0, \sigma_y^2 )$$.


$subhead p( y_t | u , theta )$$
For $latex t = 0 , \ldots , T - 1$$,
$latex y_t \sim \B{N}( u_t , \sigma_y^2 )$$.

$subhead p( u | theta )$$
For $latex t = 0$$, $latex u_t \sim \B{N}( 0 , \theta_0^2 )$$,
and for $latex t = 1 , \ldots , T - 1$$,
$latex u_t - u_{t-1} \sim \B{N}( 0 , \theta_0^2 )$$.


$head Input$$

$subhead trace_optimize_fixed$$
This is either $code true$$ or $code false$$.
If it is true, a $icode%print_level% = 5%$$
$cref/trace/ipopt_trace/$$ of the fixed effects optimization
is included in the program output.
Otherwise the ipopt $icode print_level$$ is zero and
no such trace is printed.

$subhead ipopt_solve$$
This is either $code true$$ or $code false$$.
If it is true, the $code CppAD::ipopt::solve$$
routine is used for optimizing the random effects,
otherwise $code CppAD::mixed::ipopt_random$$ is used; see
$cref/evaluation_method/optimize_random/options/evaluation_method/$$.

$subhead random_seed$$
This is a non-negative integer equal to the
seed for the random number generator,
to be specific,
$cref/s_in/manage_gsl_rng/new_gsl_rng/s_in/$$ used during the call to
$code new_gsl_rng$$.

$subhead number_times$$
This is a positive integer specifying the number of time points.
This is also the number of random effects and number of data values.

$head Output$$

$subhead actual_seed$$
If $icode random_seed$$ is zero,
the system clock, instead of $icode random_seed$$,
is used to seed the random number generator.
The actual random seed $icode actual_seed$$ is printed
so that you can reproduce results when $icode random_seed$$ is zero.

$subhead initialize_bytes$$
Is the amount of heap memory, in bytes,
added to the derived class object
during its $cref initialize$$ call.
Note that more temporary memory may have been used during this call.

$subhead initialize_seconds$$
Is the number of seconds used by the derived class $cref initialize$$ call.

$subhead optimize_fixed_seconds$$
Is the number of seconds used by the call to
$cref optimize_fixed$$ that is used to compute the
optimal fixed effects.

$subhead optimize_random_seconds$$
Is the number of seconds used by a single call to
$cref optimize_random$$ that is used to compute the
optimal random effects.

$subhead information_mat_seconds$$
Is the number of seconds used by the call to
$cref information_mat$$ that computes the observed information matrix.

$subhead sample_fixed_seconds$$
Is the number of seconds used by the call to
$cref sample_fixed$$ that computes the
$cref/number_sample_fixed/capture_xam.cpp/Input/number_fixed_samples/$$
samples for the fixed effects.

$subhead theta_0_estimate$$
Is the optimal estimate for $latex \theta_0$$; see the
$cref/problem/ar1_xam.cpp/Problem/$$ definition.

$subhead ar1_xam_ok$$
If this program passes it's correctness test,
$icode ar1_xam_ok$$ is true and the program return code is $code 0$$.
Otherwise $icode ar1_xam_ok$$ it is false and the return code is $code 1$$.

$head Example$$
The $code cppad_mixed$$ automated testing system uses the following
values for the arguments to $code ar1_xam$$:
$code
$verbatim%speed/CMakeLists.txt
	%0%# BEGIN ar1_xam arguments%# END ar1_xam arguments%0%$$
$$

$head Source Code$$
$code
$srcfile%speed/ar1_xam.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$


$end
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <gsl/gsl_randist.h>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/sparse_mat_info.hpp>
# include <cppad/mixed/manage_gsl_rng.hpp>
# include <Eigen/Dense>

namespace {
	using std::cout;
	using std::string;
	using CppAD::vector;
	using CppAD::log;
	using CppAD::AD;
	using CppAD::mixed::sparse_mat_info;
	//
	typedef AD<double>    a1_double;
	typedef AD<a1_double> a2_double;
	//
	// The constant sigma_y
	const double sigma_y = 0.1;
	//
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
		{	assert( n_fixed == 1);
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
			//
			// initialize summation
			vec[0] = a2_double(0.0);
			//
			// sqrt_2pi = CppAD::sqrt(8.0 * CppAD::atan(1.0) );
			//
			for(size_t i = 0; i < n_random_; i++)
			{	a2_double sigma  = a2_double(sigma_y);
				a2_double res    = (y_[i] - u[i]) / sigma;
				//
				// p(y_i | u, theta)
				vec[0] += res * res / a2_double(2.0);
				// following term does not depend on fixed or random effects
				// vec[0] += log(sqrt_2pi);
				//
				// p(u_i | theta)
				sigma = theta[0];
				if( i == 0 )
					res = u[i] / sigma;
				else
					res = ( u[i] - u[i-1] ) / sigma;
				vec[0] += log(sigma) + res * res / a2_double(2.0);
				// following term does not depend on fixed or random effects
				// vec[0] += log(sqrt_2pi);
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
		size_t n_digits = 5 + size_t( std::log10(value) + 1e-9 );
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
		"ipopt_solve",
		"random_seed",
		"number_times"
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
	assert( n_arg == 4 );
	bool   trace_optimize_fixed   = string( argv[1] ) == "true";
	bool   ipopt_solve            = string( argv[2] ) == "true";
	size_t random_seed            = std::atoi( argv[3] );
	size_t number_times           = std::atoi( argv[4] );
	//
	// print the command line arugments with labels for each value
	for(size_t i = 0; i < n_arg; i++)
		label_print(arg_name[i], argv[1+i]);
	//
	// initialize gsl random number generator
	size_t actual_seed = CppAD::mixed::new_gsl_rng(random_seed);
	label_print("actual_seed", actual_seed);
	//
	size_t n_data   = number_times;
	size_t n_random = number_times;
	size_t n_fixed  = 1;
	vector<double>
		fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
	fixed_lower[0] = 1e-5; fixed_in[0] = 2.0; fixed_upper[0] = inf;
	//
	// explicit constriants
	vector<double> fix_constraint_lower(0), fix_constraint_upper(0);
	//
	//
	// difference
	gsl_rng* rng = CppAD::mixed::get_gsl_rng();
	vector<double> y(n_data), random_in(n_random), theta_sim(n_fixed);
	theta_sim[0] = 1.0;
	for(size_t i = 0; i < n_data; i++)
	{	y[i]         = double(i + 1) * theta_sim[0];
		y[i]        += gsl_ran_gaussian(rng, sigma_y);
		random_in[i] = 0.0;
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed = true;
	CppAD::mixed::sparse_mat_info A_info; // empty matrix
	mixed_derived mixed_object(n_fixed, n_random, quasi_fixed, A_info, y);

	// initialization
	double start_seconds = CppAD::elapsed_seconds();
	std::map<string, size_t> size_map =
		mixed_object.initialize(fixed_in, random_in);
	double end_seconds = CppAD::elapsed_seconds();
	//
	// print amoumt of memory added to mixed_object during initialize
	// (use commans to separate every three digits).
	size_t num_bytes_before = size_map["num_bytes_before"];
	size_t num_bytes_after  = size_map["num_bytes_after"];
	string bytes_str = CppAD::to_string(num_bytes_after - num_bytes_before);
	size_t n_char = bytes_str.size();
	string initialize_bytes = "";
	for(size_t i = 0; i < n_char; i++)
	{	if( i != 0 && (n_char - i) % 3 == 0 )
			initialize_bytes += ",";
		initialize_bytes += bytes_str[i];
	}
	label_print("initialize_bytes", initialize_bytes);
	label_print("initialize_seconds", end_seconds - start_seconds);

	// optimize the fixed effects using quasi-Newton method
	string random_ipopt_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           none\n"
		"Numeric tol                       1e-7\n"
	;
	if( ipopt_solve ) random_ipopt_options +=
		"String  evaluation_method         ipopt_solve\n";
	string fixed_ipopt_options =
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
	if( trace_optimize_fixed )
		std::cout << std::endl;
	label_print("optimize_fixed_seconds", end_seconds - start_seconds);
	//
	// estimate of fixed effects
	vector<double> theta_out = solution.fixed_opt;
	//
	// estimate of random effects
	start_seconds = CppAD::elapsed_seconds();
	vector<double> u_out     = mixed_object.optimize_random(
		random_ipopt_options,
		theta_out,
		random_lower,
		random_upper,
		random_in
	);
	end_seconds = CppAD::elapsed_seconds();
	label_print("optimize_random_seconds", end_seconds - start_seconds);
	//
	// information matrix
	start_seconds = CppAD::elapsed_seconds();
	CppAD::mixed::sparse_mat_info
	information_info = mixed_object.information_mat(solution, u_out);
	end_seconds = CppAD::elapsed_seconds();
	label_print("information_mat_seconds", end_seconds - start_seconds);
	//
	// sample approximate posteroior for fixed effects
	size_t number_fixed_samples = 1000;
	vector<double> sample( number_fixed_samples * n_fixed );
	start_seconds = CppAD::elapsed_seconds();
	mixed_object.sample_fixed(
		sample,
		information_info,
		solution,
		fixed_lower,
		fixed_upper,
		u_out
	);
	end_seconds = CppAD::elapsed_seconds();
	label_print("sample_fixed_seconds", end_seconds - start_seconds);
	//
	// compute the sample standard deviations
	// and ratio of error divided by sample standard deviation
	vector<double> sample_std(n_fixed), estimate_ratio(n_fixed);
	for(size_t j = 0; j < n_fixed; j++)
		sample_std[j] = 0.0;
	for(size_t i = 0; i < number_fixed_samples; i++)
	{	for(size_t j = 0; j < n_fixed; j++)
		{	double diff    = sample[i * n_fixed + j] - theta_out[j];
			sample_std[j] += diff * diff;
		}
	}
	for(size_t j = 0; j < n_fixed; j++)
	{	sample_std[j] = std::sqrt( sample_std[j] / number_fixed_samples );
		// check if results results are reasonable
		estimate_ratio[j] = ( theta_out[j] - theta_sim[j] ) / sample_std[j];
		ok  &= std::fabs(estimate_ratio[j]) < 4.0;
	}
	label_print("theta_0_estimate", theta_out[0] );
	label_print("theta_0_std",      sample_std[0] );
	label_print("theta_0_ratio",    estimate_ratio[0] );
	//
	// Check that only the lower limit on theta[1] is active.
	ok &= solution.fixed_lag.size() == n_fixed;
	ok &= solution.fixed_lag[0] == 0.0;
	ok &= solution.fix_con_lag.size() == 0;
	ok &= solution.ran_con_lag.size() == 0;
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
