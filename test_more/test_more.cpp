// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <iostream>
# include <cassert>
# include <cstring>

// cppad_mixed subdirectory
extern bool abs_fix_con(void);
extern bool binomial(void);
extern bool delta_ran_obj(void);
extern bool der_var_hes(void);
extern bool fixed_lag(void);
extern bool ldlt_cholmod(void);
extern bool max_iter_neg(void);
extern bool n_mixture(void);
extern bool no_fix_likelihood(void);
extern bool no_random(void);
extern bool no_random_info(void);
extern bool opt_ran_fail(void);
extern bool ran_obj_tst(void);
extern bool laplace_obj_hes(void);
extern bool laplace_obj_tst(void);
extern bool sample_fixed_1(void);
extern bool sample_fixed_2(void);
extern bool scale(void);
extern bool solution_check(void);
extern bool zero_random_one(void);
extern bool zero_random_two(void);

// anonymous namespace
namespace {
	using std::cout;
	using std::endl;

	// function that runs one test
	static size_t Run_ok_count    = 0;
	static size_t Run_error_count = 0;
	void Run(bool test_fun(void), const char* test_name)
	{
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
	// This comment expected by bin/test_one.sh
	RUN(abs_fix_con);
	RUN(binomial);
	RUN(delta_ran_obj);
	RUN(der_var_hes);
	RUN(fixed_lag);
	RUN(ldlt_cholmod);
	RUN(max_iter_neg);
	RUN(n_mixture);
	RUN(no_fix_likelihood);
	RUN(no_random);
	RUN(no_random_info);
	RUN(opt_ran_fail);
	RUN(ran_obj_tst);
	RUN(laplace_obj_hes);
	RUN(laplace_obj_tst);
	RUN(sample_fixed_1);
	RUN(sample_fixed_2);
	RUN(scale);
	RUN(solution_check);
	RUN(zero_random_one);
	RUN(zero_random_two);
	// This comment also expected by bin/test_one.sh

	// summary report
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
