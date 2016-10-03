// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin manage_gsl_rng.cpp$$
$spell
	CppAD
	gsl
	rng
$$

$section Manage GSL Random Number Generator: Example and Test$$

$code
$srcfile%example/user/manage_gsl_rng.cpp%
	0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/utility.hpp> // CppAD::vector
# include <gsl/gsl_randist.h>
# include <cppad/mixed/manage_gsl_rng.hpp>

bool manage_gsl_rng_xam(void)
{	bool ok = true;
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
	{	ok    &= 0. <= sim[i];
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
	{	for(j = 0; j < i; j++)
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
