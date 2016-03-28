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

extern bool abs_density_xam(void);
extern bool cholmod_solve2_xam(void);
extern bool cholmod_xam(void);
extern bool data_mismatch_xam(void);
extern bool derived_xam(void);
extern bool eigen_xam(void);
extern bool fix_con_eval_xam(void);
extern bool fix_con_hes_xam(void);
extern bool fix_con_jac_xam(void);
extern bool fix_constraint_xam(void);
extern bool fix_like_eval_xam(void);
extern bool fix_like_hes_xam(void);
extern bool fix_like_jac_xam(void);
extern bool hes_cross_xam(void);
extern bool ipopt_xam_run(void);
extern bool lasso_xam(void);
extern bool logdet_jac_xam(void);
extern bool manage_gsl_rng_xam(void);
extern bool newton_step_xam(void);
extern bool no_random_xam(void);
extern bool optimize_fixed_xam(void);
extern bool optimize_random_xam(void);
extern bool ran_con_eval_xam(void);
extern bool ran_con_jac_xam(void);
extern bool ran_constraint_xam(void);
extern bool ran_hes_fun_xam(void);
extern bool ran_like_jac_xam(void);
extern bool ran_likelihood_xam(void);
extern bool ran_objcon_hes_xam(void);
extern bool ran_obj_eval_xam(void);
extern bool ran_obj_jac_xam(void);
extern bool update_factor_xam(void);

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
	RUN(abs_density_xam);
	RUN(cholmod_solve2_xam);
	RUN(cholmod_xam);
	RUN(data_mismatch_xam);
	RUN(derived_xam);
	RUN(eigen_xam);
	RUN(fix_con_eval_xam);
	RUN(fix_con_hes_xam);
	RUN(fix_con_jac_xam);
	RUN(fix_constraint_xam);
	RUN(fix_like_eval_xam);
	RUN(fix_like_hes_xam);
	RUN(fix_like_jac_xam);
	RUN(hes_cross_xam);
	RUN(ipopt_xam_run);
	RUN(lasso_xam);
	RUN(logdet_jac_xam);
	RUN(manage_gsl_rng_xam);
	RUN(newton_step_xam);
	RUN(no_random_xam);
	RUN(optimize_fixed_xam);
	RUN(optimize_random_xam);
	RUN(ran_con_eval_xam);
	RUN(ran_con_jac_xam);
	RUN(ran_constraint_xam);
	RUN(ran_hes_fun_xam);
	RUN(ran_like_jac_xam);
	RUN(ran_likelihood_xam);
	RUN(ran_objcon_hes_xam);
	RUN(ran_obj_eval_xam);
	RUN(ran_obj_jac_xam);
	RUN(update_factor_xam);

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
