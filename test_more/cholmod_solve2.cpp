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
$begin cholmod_solve2_test$$
$spell
	cholmod
$$

$section Test Where cholmod_solve2 Does not Detect a Non-positive Matrix$$

$head matrix$$
$latex \[
	A = \left( \begin{array}{cc}
		1 & 0 \\
		0 & -1
	\end{array} \right)
\] $$

$head Problem$$
The following check in the example below returns $icode%ok% = true%$$:
$srccode%cpp%
	// One might expect com.status to be equal to CHOLMOD_NOT_POSDEF
	ok &= ! (com.status == CHOLMOD_NOT_POSDEF);
%$$

$head Source Code$$
$srcfile%test_more/cholmod_solve2.cpp%0%// BEGIN C++%// END C++%1%$$

$end
*/
// BEGIN C++
# include <cholmod.h>
# include <limits>
# include <cmath>
# include <cassert>

# define CHOLMOD_TRUE                  1
# define CHOLMOD_FALSE                 0
# define CHOLMOD_STYPE_LOWER_TRIANGLE -1

# define DEMONSTRATE_PROBLEM 0

namespace { // BEGIN_EMPTY_NAMESPACE
	void add_T_entry(cholmod_triplet *T, int r, int c, double x)
	{	size_t k           = T->nnz;
		((int*)T->i)[k]    = r;
		((int*)T->j)[k]    = c;
		((double*)T->x)[k] = x;
		T->nnz++;
	}
}

bool cholmod_solve2(void)
{	bool ok = true;

	// Elements for the the matrix A
	size_t nrow = 2;
	double A_data[] = {
		 1.0,   0.0,
		 0.0,  -1.0
	};

	// initialize common
	cholmod_common com;
	cholmod_start(&com);

	// always do simplicial factorization
	com.supernodal = CHOLMOD_SIMPLICIAL;

	// do LDL' factorization and leave in LDL' form
	com.final_ll = CHOLMOD_FALSE;

	// allocate triplet
	size_t ncol = nrow;
	size_t nzmax = nrow * ncol;
	int    T_stype = CHOLMOD_STYPE_LOWER_TRIANGLE;
	int    T_xtype = CHOLMOD_REAL;
	cholmod_triplet *T =
		cholmod_allocate_triplet(nrow, ncol, nzmax, T_stype, T_xtype, &com);
	ok &= T->nnz ==  0;

	// triplet entries corresponding to lower triangle of A
	for(size_t i = 0; i < nrow; i++)
	{	for(size_t j = 0; j <= i; j++)
			add_T_entry(T, i, j, A_data[i * ncol + j] );
	}
	ok &= T->nnz == 3;

	// convert triplet to sparse representation of A
	cholmod_sparse* A = cholmod_triplet_to_sparse(T, 0, &com);

	// factor the matrix
	cholmod_factor *L = cholmod_analyze(A, &com);
	int flag = cholmod_factorize(A, L, &com);

	// check properties of factor
	assert( flag     == CHOLMOD_TRUE );  // return flag OK
	assert( L->n     == nrow );          // number of rows and coluns
	assert( L->minor == nrow );          // successful factorization
	assert( L->is_ll == CHOLMOD_FALSE ); // factorization is LDL'

	// ----------------------------------------------------------------------
	// One might expect com.status to be equal to CHOLMOD_NOT_POSDEF
	ok &= ! (com.status == CHOLMOD_NOT_POSDEF);
	// ----------------------------------------------------------------------

	// free memory
	cholmod_free_triplet(&T,    &com);
	cholmod_free_sparse( &A,    &com);
	cholmod_free_factor( &L,    &com);

	// finish up
	cholmod_finish(&com);

	// check count of malloc'ed - free'd
	ok &= com.malloc_count == 0;

	return ok;
}
// END C++
