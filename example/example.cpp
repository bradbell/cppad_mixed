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
$begin example_devel.cpp$$

$sectin Run C++ Examples$$

$head Syntax$$
$code example/devel/cpp$$

$head Purpose$$
This runs all the C++ examples and prints out their correctness
test results together with a summary result at the end.

$end
-----------------------------------------------------------------------------
*/
# include <iostream>
# include <cassert>
# include <cstring>
# include <cppad/mixed/configure.hpp>

// optional
extern bool cholmod_xam(void);

// cppad_mixed subdirectory
extern bool fix_constraint_xam(void);
extern bool derived_xam(void);
extern bool fix_constraint_eval_xam(void);
extern bool fix_constraint_hes_xam(void);
extern bool fix_constraint_jac_xam(void);
extern bool data_mismatch_xam(void);
extern bool logdet_grad_xam(void);
extern bool ranobj_grad_xam(void);
extern bool eigen_xam(void);
extern bool fix_like_hes_xam(void);
extern bool fix_like_jac_xam(void);
extern bool hes_cross_xam(void);
extern bool hes_ran_fun_xam(void);
extern bool ranobj_eval_xam(void);
extern bool ipopt_xam_run(void);
extern bool newton_step_xam(void);
extern bool no_random_xam(void);
extern bool optimize_fixed_xam(void);
extern bool optimize_random_xam(void);
extern bool fix_like_eval_xam(void);
extern bool ran_likelihood_grad_xam(void);
extern bool ranobj_hes_xam(void);
extern bool manage_gsl_rng_xam(void);

// anonymous namespace
namespace {
	// function that runs one test
	static size_t Run_ok_count    = 0;
	static size_t Run_error_count = 0;
	void Run(bool test_fun(void), const char* test_name)
	{	using std::cout;
		using std::endl;
		//
		std::streamsize width = 30;
		cout.width( width );
		cout.setf( std::ios_base::left );
		cout << test_name << ':';
		assert( std::strlen(test_name) < size_t(width) );
		//
		bool ok = test_fun();
		if( ok )
		{	cout << "OK" << endl;
			Run_ok_count++;
		}
		else
		{	cout << "Error" << endl;
			Run_error_count++;
		}
	}
}
// macro for calls Run
# define RUN(test_name) Run( test_name, #test_name )

// main program that runs all the tests
int main(void)
{
	RUN(fix_constraint_xam);
	RUN(derived_xam);
	RUN(fix_constraint_eval_xam);
	RUN(fix_constraint_hes_xam);
	RUN(fix_constraint_jac_xam);
	RUN(data_mismatch_xam);
	RUN(logdet_grad_xam);
	RUN(ranobj_grad_xam);
	RUN(ranobj_hes_xam);
	RUN(eigen_xam);
	RUN(fix_like_hes_xam);
	RUN(fix_like_jac_xam);
	RUN(hes_cross_xam);
	RUN(hes_ran_fun_xam);
	RUN(ranobj_eval_xam);
	RUN(ipopt_xam_run);
	RUN(manage_gsl_rng_xam);
	RUN(newton_step_xam);
	RUN(ran_likelihood_grad_xam);
	RUN(no_random_xam);
	RUN(optimize_fixed_xam);
	RUN(optimize_random_xam);
	RUN(fix_like_eval_xam);


# if CPPAD_MIXED_HAS_SUITESPARSE
	RUN(cholmod_xam);
# endif

	// summary report
	using std::cout;
	using std::endl;
	int return_flag;
	if( Run_error_count == 0 )
	{	cout << "All " << Run_ok_count << " tests passed." << endl;
		return_flag = 0;
	}
	else
	{	cout << Run_error_count << " tests failed." << endl;
		return_flag = 1;
	}

	return return_flag;
}
