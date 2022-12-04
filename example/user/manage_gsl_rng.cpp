// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin manage_gsl_rng.cpp}

Manage GSL Random Number Generator: Example and Test
####################################################

{xrst_literal
   // BEGIN C++
   // END C++
}

{xrst_end manage_gsl_rng.cpp}
*/
// BEGIN C++
# include <cppad/utility/vector.hpp>
# include <gsl/gsl_randist.h>
# include <cppad/mixed/manage_gsl_rng.hpp>

bool manage_gsl_rng_xam(void)
{  bool ok = true;
   size_t s_in = 123;
   size_t s_out;

   // initialize random number generator using a specific seed
   s_out = CppAD::mixed::new_gsl_rng(s_in);
   ok   &= s_out == s_in;
   size_t i, j, n = 10;
   CppAD::vector<double> sim(n);
   for(i = 0; i < n; i++)
      sim[i] = gsl_ran_flat(CppAD::mixed::get_gsl_rng(), 0., 1.);
   for(i = 0; i < n; i++)
   {  ok    &= 0. <= sim[i];
      ok    &= sim[i] <= 1.;
      for(j = 0; j < n; j++)
         ok &= (i == j) || (sim[i] != sim[j]);
   }
   // done with this random number generator
   CppAD::mixed::free_gsl_rng();

   // test running the same seed
   CppAD::mixed::new_gsl_rng(s_out);
   for(i = 0; i < n; i++)
      ok    &= ( sim[i] == gsl_ran_flat(CppAD::mixed::get_gsl_rng(), 0., 1.) );

   // done with this random number generator
   CppAD::mixed::free_gsl_rng();

   // test using system time for the seed
   s_in  = 0;
   s_out = CppAD::mixed::new_gsl_rng(s_in);
   ok   &= s_out != s_in;
   CppAD::vector<double> temp(n);
   for(i = 0; i < n; i++)
      temp[i] = gsl_ran_flat(CppAD::mixed::get_gsl_rng(), 0., 1.);
   for(i = 0; i < n; i++)
   {  for(j = 0; j < i; j++)
         ok &= (temp[i] != sim[j]);
   }

   // make this the previous simulation
   sim = temp;

   // done with this random number generator
   CppAD::mixed::free_gsl_rng();

   // test using a different system time for the seed
   CppAD::mixed::new_gsl_rng(0);
   for(i = 0; i < n; i++)
      ok &= ( sim[i] != gsl_ran_flat(CppAD::mixed::get_gsl_rng(), 0., 1.) );

   // done with this random number generator
   CppAD::mixed::free_gsl_rng();

   // test using the random seed chosen automatically
   CppAD::mixed::new_gsl_rng(s_out);
   for(i = 0; i < n; i++)
      ok &= ( sim[i] == gsl_ran_flat(CppAD::mixed::get_gsl_rng(), 0., 1.) );

   // done with this random number generator
   CppAD::mixed::free_gsl_rng();

   return ok;
}
// END C++
