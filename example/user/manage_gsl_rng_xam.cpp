// $Id:$
/* --------------------------------------------------------------------------
dismod_at: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin manage_gsl_rng_xam.cpp$$
$spell
	gsl
	rng
$$

$section Manage GSL Random Number Generator: Example and Test$$

$code
$verbatim%example/devel/utility/manage_gsl_rng_xam.cpp%
	0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/vector.hpp>
# include <gsl/gsl_randist.h>
# include <dismod_at/manage_gsl_rng.hpp>

bool manage_gsl_rng_xam(void)
{	bool ok = true;
	size_t s_in = 123;
	size_t s_out;

	// initialize random number generator using a specific seed
	s_out = dismod_at::new_gsl_rng(s_in);
	ok   &= s_out == s_in;
	size_t i, j, n = 10;
	CppAD::vector<double> sim(n);
	for(i = 0; i < n; i++)
		sim[i] = gsl_ran_flat(dismod_at::get_gsl_rng(), 0., 1.);
	for(i = 0; i < n; i++)
	{	ok    &= 0. <= sim[i];
		ok    &= sim[i] <= 1.;
		for(j = 0; j < n; j++)
			ok &= (i == j) || (sim[i] != sim[j]);
	}
	// done with this random number generator
	dismod_at::free_gsl_rng();

	// test running the same seed
	dismod_at::new_gsl_rng(s_out);
	for(i = 0; i < n; i++)
		ok    &= ( sim[i] == gsl_ran_flat(dismod_at::get_gsl_rng(), 0., 1.) );

	// done with this random number generator
	dismod_at::free_gsl_rng();

	// test using system time for the seed
	s_in  = 0;
	s_out = dismod_at::new_gsl_rng(s_in);
	ok   &= s_out != s_in;
	CppAD::vector<double> temp(n);
	for(i = 0; i < n; i++)
		temp[i] = gsl_ran_flat(dismod_at::get_gsl_rng(), 0., 1.);
	for(i = 0; i < n; i++)
	{	for(j = 0; j < i; j++)
			ok &= (temp[i] != sim[j]);
	}

	// make this the previous simulation
	sim = temp;

	// done with this random number generator
	dismod_at::free_gsl_rng();

	// test using a different system time for the seed
	dismod_at::new_gsl_rng(0);
	for(i = 0; i < n; i++)
		ok &= ( sim[i] != gsl_ran_flat(dismod_at::get_gsl_rng(), 0., 1.) );

	// done with this random number generator
	dismod_at::free_gsl_rng();

	// test using the random seed chosen automatically
	dismod_at::new_gsl_rng(s_out);
	for(i = 0; i < n; i++)
		ok &= ( sim[i] == gsl_ran_flat(dismod_at::get_gsl_rng(), 0., 1.) );

	// done with this random number generator
	dismod_at::free_gsl_rng();

	return ok;
}
// END C++
