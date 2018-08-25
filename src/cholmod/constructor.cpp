// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin ldlt_cholmod_ctor$$
$spell
	ldlt
	nrow
	xam
	ldlt_obj
	cholmod
	Cpp
	chol
	hes
	initializes
	simplicial
	supernodal
$$

$section Cholmod LDLT Constructor$$

$head Syntax$$
$codei%CppAD::mixed::ldlt_cholmod %ldlt_obj%(%nrow%)%$$


$head Private$$
The $cref ldlt_cholmod$$ class is an
$cref/implementation detail/ldlt_cholmod/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head nrow_$$
The argument $icode nrow$$ has prototype
$codei%
	size_t %nrow%
%$$
It is the number of rows in the symmetric matrix
we will compute the LDLT factor of.
The member variable $code nrow_$$ is set to this value.

$head Pointers$$
This sets all the pointer private member variable equal to null.

$head common_$$
The $code common_$$ variable is initialized as follows:
$codei%
	cholmod_start(&common_);
%$$

$subhead Simplicial Factorization$$
The following value is set (and does not change):
$codei%
	common_.supernodal = CHOLMOD_SIMPLICIAL;
%$$

$subhead LDL' Factorization$$
The following value is set (and does not change):
$codei%
	common_.final_ll = CHOLMOD_FLASE;
%$$

$head Example$$
The file $cref/ldlt_cholmod.cpp/ldlt_cholmod.cpp/constructor/$$ contains an
example and test that uses this constructor.

$end
*/

# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cassert>
# include <cppad/cppad.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

ldlt_cholmod::ldlt_cholmod(size_t nrow)
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
{	assert( CPPAD_NULL == NULL );
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
$begin ldlt_cholmod_dtor$$
$spell
	ldlt
	Cholmod
	Cpp
$$

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
ldlt_cholmod::~ldlt_cholmod(void)
{	// free all the private pointers
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
