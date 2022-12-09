// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin ldlt_cholmod_ctor dev}
{xrst_spell
   ll
   nrow
   simplicial
   supernodal
}

Cholmod LDLT Constructor
########################

Syntax
******
``CppAD::mixed::ldlt_cholmod`` *ldlt_obj* ( *nrow* )

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE
   // END_PROTOTYPE
}

Private
*******
The :ref:`ldlt_cholmod-name` class is an
:ref:`implementation detail<ldlt_cholmod@Private>` and not part of the
CppAD Mixed user API.

nrow\_
******
The argument *nrow*
is the number of rows in the symmetric matrix
we will compute the LDLT factor of.
The member variable ``nrow_`` is set to this value.

Pointers
********
This sets all the pointer private member variable equal to null.

common\_
********
The ``common_`` variable is initialized as follows:

   ``cholmod_start`` (& ``common_`` );

Simplicial Factorization
========================
The following value is set (and does not change):

   ``common_.supernodal`` = ``CHOLMOD_SIMPLICIAL`` ;

LDL' Factorization
==================
The following value is set (and does not change):

   ``common_.final_ll`` = ``CHOLMOD_FLASE`` ;

Example
*******
The file :ref:`ldlt_cholmod.cpp<ldlt_cholmod.cpp@constructor>` contains an
example and test that uses this constructor.

{xrst_end ldlt_cholmod_ctor}
*/

# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cassert>
# include <cppad/cppad.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN_PROTOTYPE
ldlt_cholmod::ldlt_cholmod(size_t nrow)
// END_PROTOTYPE
:
nrow_          (nrow)            ,
init_done_     (false)           ,
update_called_ (false)           ,
sym_matrix_    (CPPAD_NULL)      ,
factor_        (CPPAD_NULL)      ,
rhs_           (CPPAD_NULL)      ,
rhs_set_       (CPPAD_NULL)      ,
sol_           (CPPAD_NULL)      ,
sol_set_       (CPPAD_NULL)      ,
work_one_      (CPPAD_NULL)      ,
work_two_      (CPPAD_NULL)
{  assert( CPPAD_NULL == NULL );
   cholmod_start(&common_);

   // Both simplical and supernodal have been tested.
   // common_.supernodal = CHOLMOD_SIMPLICIAL;
   // common_.supernodal = CHOLMOD_AUTO;
   // common_.supernodal = CHOLMOD_SUPERNODAL;
   common_.supernodal = CHOLMOD_SIMPLICIAL;

   // is factorization LLT or LDLT
   common_.final_ll = CHOLMOD_FALSE;

}

/*
{xrst_begin ldlt_cholmod_dtor dev}

Cholmod Destructor
##################

Private
*******
This class is an implementation detail and not part of the
CppAD Mixed user API.

Discussion
**********
This frees the memory corresponding to all the
private member variables and the clear the common with

   ``cholmod_finish`` (& ``common_`` );

{xrst_end ldlt_cholmod_dtor}
*/
ldlt_cholmod::~ldlt_cholmod(void)
{  // free all the private pointers
   cholmod_free_sparse (&sym_matrix_, &common_ );
   cholmod_free_factor (&factor_,     &common_ );
   cholmod_free_dense  (&rhs_,        &common_ );
   cholmod_free_sparse (&rhs_set_,    &common_ );
   cholmod_free_dense  (&sol_,        &common_ );
   cholmod_free_sparse (&sol_set_,    &common_ );
   cholmod_free_dense  (&work_one_,   &common_ );
   cholmod_free_dense  (&work_two_,   &common_ );
   // clear common
   cholmod_finish(&common_);

   // always do simplicial factorization
   assert( common_.supernodal == CHOLMOD_SIMPLICIAL );

   // do LDL' factorization and leave in LDL' form
   assert( common_.final_ll == CHOLMOD_FALSE );

   // check nothing left allocated
   assert( common_.malloc_count == 0 );
}


} } // END_CPPAD_MIXED_NAMESPACE
