// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */

# include <cppad/mixed/cholmod.hpp>
# include <cassert>
# include <cppad/cppad.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

/*
$begin cholmod_ctor$$

$section Cholmod Constructor$$

$head Syntax$$
$codei%CppAD::mixed::cholmod %chol_ran_hes%;$$

$head Private$$
This class is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Discussion$$
This sets all the private member variables equal to the null pointer
and initializes the common block with the call
$codei%
	cholmod_start(&common_);
%$$

$end
*/
cholmod::cholmod(void)
:
triplet_   (CPPAD_NULL)      ,
pos_matrix_(CPPAD_NULL)      ,
factor_    (CPPAD_NULL)      ,
rhs_       (CPPAD_NULL)      ,
rhs_set_   (CPPAD_NULL)      ,
sol_       (CPPAD_NULL)      ,
sol_set_   (CPPAD_NULL)      ,
work_one_  (CPPAD_NULL)      ,
work_two_  (CPPAD_NULL)
{	assert( CPPAD_NULL == NULL );
	cholmod_start(&common_);
}

/*
$begin cholmod_dtor$$

$section Cholmod Destructor$$

$head Private$$
This class is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Discussion$$
This frees the memory corresponding to all the
private member variables and the clear the common with
$codei%
	cholmod_finish(&common_);
%$$

$end
*/
cholmod::~cholmod(void)
{	// free all the private pointers
	cholmod_free_triplet(&triplet_,    &common_ );
	cholmod_free_sparse (&pos_matrix_, &common_ );
	cholmod_free_factor (&factor_,     &common_ );
	cholmod_free_dense  (&rhs_,        &common_ );
	cholmod_free_sparse (&rhs_set_,    &common_ );
	cholmod_free_dense  (&sol_,        &common_ );
	cholmod_free_sparse (&sol_set_,    &common_ );
	cholmod_free_dense  (&work_one_,   &common_ );
	cholmod_free_dense  (&work_two_,   &common_ );
	// clear common
	cholmod_finish(&common_);

	// check nothing left allocated
	assert( common_.malloc_count == 0 );
}


} } // END_CPPAD_MIXED_NAMESPACE
