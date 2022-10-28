// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
$begin cholesky_devel.cpp$$
$spell
   cholesky
   devel
   cpp
$$

$section Run C++ Examples$$

$head Syntax$$
$code cholesky/cholesky$$

$head Purpose$$
This runs all the C++ sparse Cholesky examples (and tests)
and prints out their correctness
test results together with a summary result at the end.

$end
-----------------------------------------------------------------------------
*/
# include <iostream>
# include <cassert>
# include <cstring>

extern bool sparse_ad_chol_eval(void);
extern bool sparse_ad_chol_eq(void);
extern bool sparse_ad_chol_perm(void);
extern bool sparse_ad_chol_sp1(void);
extern bool sparse_ad_chol_sp2(void);
extern bool sparse_ad_chol_var(void);

// anonymous namespace
namespace {
   // function that runs one test
   static size_t Run_ok_count    = 0;
   static size_t Run_error_count = 0;
   void Run(bool test_fun(void), const char* test_name)
   {  using std::cout;
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
      {  cout << "OK" << endl;
         Run_ok_count++;
      }
      else
      {  cout << "Error" << endl;
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
   RUN(sparse_ad_chol_eval);
   RUN(sparse_ad_chol_eq);
   RUN(sparse_ad_chol_perm);
   RUN(sparse_ad_chol_sp1);
   RUN(sparse_ad_chol_sp2);
   RUN(sparse_ad_chol_var);
   // This comment also expected by bin/test_one.sh

   // summary report
   using std::cout;
   using std::endl;
   int return_flag;
   if( Run_error_count == 0 )
   {  cout << "All " << Run_ok_count << " tests passed." << endl;
      return_flag = 0;
   }
   else
   {  cout << Run_error_count << " tests failed." << endl;
      return_flag = 1;
   }

   return return_flag;
}
