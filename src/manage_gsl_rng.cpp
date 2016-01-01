// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin manage_gsl_rng$$
$spell
	CppAD
	std
	gsl
	rng
$$

$section Set, Get, And Free A GSL Random Number Generator$$

$head Syntax$$
$icode%s_out% = CppAD::mixed::new_gsl_rng(%s_in%)
%$$
$icode%rng% = CppAD::mixed::get_gsl_rng()
%$$
$codei%CppAD::mixed::free_gsl_rng()
%$$

$head Purpose$$
Create and use a GSL random number generator.

$head new_gsl_rng$$
This routine creates a new GSL random number generator.
If a previous random number generator was created, it must
be freed using $code free_gsl_rng$$ before $code new_gsl_rng$$
can be called again.

$subhead s_in$$
This argument has prototype
$codei%
	size_t %s_in%
%$$
If $icode%s_in% != 0%$$,
it is used as a seed for the random number generator.
Otherwise the actual seed is the number of seconds returned by
$code std::time$$ plus the number of previous calls to $code set_gsl_rng$$.
(Adding the number of calls prevents the same
seed from being used by calls that are close together in time.)

$subhead s_out$$
This return value prototype
$codei%
	size_t %s_out%
%$$
and is the actual seed that was used to initialize the random number generator.

$head get_gsl_rng$$
If we are between a call to
$code new_gsl_rng$$ and $code free_gsl_rng$$,
this routine retrieves a pointer to the current
GSL random number generator.
Otherwise it returns the null pointer.

$subhead rng$$
The return value $icode rng$$ has prototype
$codei%
	gsl_rng* %rng%
%$$

$head free_gsl_rng$$
Once you are done with a generator created by $code new_gsl_rng$$,
you should free the corresponding memory using
$codei%
	gsl_rng_free()
%$$.


$children%
	example/user/manage_gsl_rng_xam.cpp
%$$
$head Example$$
The file $cref manage_gsl_rng_xam.cpp$$ contains an example and test of
$code manage_gsl_rng$$.  It returns $code true$$, if the test passes,
and $code false$$ otherwise.

$end
-----------------------------------------------------------------------------
*/
# include <gsl/gsl_rng.h>
# include <ctime>
# include <cassert>
# include <cppad/mixed/manage_gsl_rng.hpp>

namespace {
	static size_t count_seed_       = 0;
	static gsl_rng* const null_ptr_ = 0;
	static gsl_rng*       rng_      = null_ptr_;
}
namespace CppAD { namespace mixed {
	size_t new_gsl_rng(size_t s_in)
	{	// check that we do not currently have a random number generator
		assert( rng_ == null_ptr_ );

		// create a new one
		rng_ = gsl_rng_alloc( gsl_rng_mt19937 );

		// initialize the return value
		size_t s_out = s_in;

		// first choice for seed
		unsigned long int lseed = static_cast<unsigned long int>(s_out);
		if( s_out == 0 )
		{	std::time_t* null_ptr(0);
			std::time_t seconds = std::time(null_ptr);
			s_out = static_cast<size_t>( seconds ) + count_seed_++;
			lseed = static_cast<unsigned long int>(s_out);
		}
		gsl_rng_set(rng_, lseed);
		//
		return s_out;
	}
	gsl_rng* get_gsl_rng(void)
	{	assert( rng_ != null_ptr_ );
		return rng_;
	}
	void free_gsl_rng(void)
	{	assert( rng_ != null_ptr_ );
		gsl_rng_free(rng_);
		rng_ = null_ptr_;
	}
} }
