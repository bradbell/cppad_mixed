// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin manage_gsl_rng}
{xrst_spell
  hpp
}

Set, Get, And Free A GSL Random Number Generator
################################################

Syntax
******

| # ``include <cppad/mixed/manage_gsl_rng.hpp>``
| *s_out* = ``CppAD::mixed::new_gsl_rng`` ( *s_in* )
| *rng* = ``CppAD::mixed::get_gsl_rng`` ()
| ``CppAD::mixed::free_gsl_rng`` ()

Purpose
*******
Create and use a GSL random number generator.

new_gsl_rng
***********
This routine creates a new GSL random number generator.
If a previous random number generator was created, it must
be freed using ``free_gsl_rng`` before ``new_gsl_rng``
can be called again.

s_in
====
This argument has prototype

   ``size_t`` *s_in*

If *s_in*  != 0 ,
it is used as a seed for the random number generator.
Otherwise the actual seed is the number of seconds returned by
``std::time`` plus the number of previous calls to ``set_gsl_rng`` .
(Adding the number of calls prevents the same
seed from being used by calls that are close together in time.)

s_out
=====
This return value prototype

   ``size_t`` *s_out*

and is the actual seed that was used to initialize the random number generator.

get_gsl_rng
***********
If we are between a call to
``new_gsl_rng`` and ``free_gsl_rng`` ,
this routine retrieves a pointer to the current
GSL random number generator.
Otherwise it returns the null pointer.

rng
===
The return value *rng* has prototype

   ``gsl_rng`` * *rng*

free_gsl_rng
************
Once you are done with a generator created by ``new_gsl_rng`` ,
you should free the corresponding memory using

   ``gsl_rng_free`` ()

.
{xrst_toc_hidden
   example/user/manage_gsl_rng.cpp
}
Example
*******
The file :ref:`manage_gsl_rng.cpp-name` contains an example and test of
``manage_gsl_rng`` .  It returns ``true`` , if the test passes,
and ``false`` otherwise.

{xrst_end manage_gsl_rng}
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
   {  // check that we do not currently have a random number generator
      assert( rng_ == null_ptr_ );

      // create a new one
      rng_ = gsl_rng_alloc( gsl_rng_mt19937 );

      // initialize the return value
      size_t s_out = s_in;

      // first choice for seed
      unsigned long int lseed = static_cast<unsigned long int>(s_out);
      if( s_out == 0 )
      {  std::time_t* null_ptr(0);
         std::time_t seconds = std::time(null_ptr);
         s_out = static_cast<size_t>( seconds ) + count_seed_++;
         lseed = static_cast<unsigned long int>(s_out);
      }
      gsl_rng_set(rng_, lseed);
      //
      return s_out;
   }
   gsl_rng* get_gsl_rng(void)
   {  assert( rng_ != null_ptr_ );
      return rng_;
   }
   void free_gsl_rng(void)
   {  assert( rng_ != null_ptr_ );
      gsl_rng_free(rng_);
      rng_ = null_ptr_;
   }
} }
