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
$begin cholmod_solve_xam$$
$spell
	cholmod
$$

$section Example Using cholmod_solve With a Non-positive Matrix$$

$head Problem Description$$
Solve for $latex x$$ in the equation $latex H x = b$$ where $latex H$$
is defined below and $latex b$$ is a column of the identity matrix.
Hence the solution $latex x$$ is the corresponding column of
$latex H^{-1}$$ which for this case is equal to $latex H$$.
$latex \[
	H = \left( \begin{array}{cc}
		1 & 0 \\
		0 & -1
	\end{array} \right)
\] $$

$head Source Code$$
$srcfile%example/private/cholmod_solve.cpp%0%// BEGIN C++%// END C++%1%$$

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


namespace { // BEGIN_EMPTY_NAMESPACE
	void add_T_entry(cholmod_triplet *T, int r, int c, double x)
	{	size_t k           = T->nnz;
		((int*)T->i)[k]    = r;
		((int*)T->j)[k]    = c;
		((double*)T->x)[k] = x;
		T->nnz++;
	}
}

bool cholmod_solve_xam(void)
{	bool ok    = true;
	double eps = 10. * std::numeric_limits<double>::epsilon();

	// Elements for the the matrix H
	size_t nrow = 2;
	double H_data[] = {
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

	// triplet entries corresponding to lower triangle of H
	for(size_t i = 0; i < nrow; i++)
	{	for(size_t j = 0; j <= i; j++)
			add_T_entry(T, i, j, H_data[i * ncol + j] );
	}
	ok &= T->nnz == 3;

	// convert triplet to sparse representation of H
	cholmod_sparse* H = cholmod_triplet_to_sparse(T, 0, &com);

	// factor the matrix
	cholmod_factor *L = cholmod_analyze(H, &com);
# ifndef NDEBUG
	int flag =
# endif
	cholmod_factorize(H, L, &com);

	// check properties of factor
	assert( flag     == CHOLMOD_TRUE );  // return flag OK
	assert( L->n     == nrow );          // number of rows and coluns
	assert( L->minor == nrow );          // successful factorization
	assert( L->is_ll == CHOLMOD_FALSE ); // factorization is LDL'
	assert( com.status == CHOLMOD_OK  ); // no problem with factorization

	// initialize right hand side vector
	cholmod_dense *B = cholmod_zeros(nrow, 1, T_xtype, &com);
	double* B_x = (double *) B->x;

	for(size_t j = 0; j < ncol; j++)
	{	// solve equation with B corresponding to j-th column of identity
		B_x[j] = 1.0;
		int sys = CHOLMOD_A;
		cholmod_dense *X = cholmod_solve(sys, L, B, &com);
		double* X_x = (double *) X->x;

		// restore B
		B_x[j] = 0.0;

		// check solution
		for(size_t i = 0; i < nrow; i++)
		{	// inv(H) = H
			double check = H_data[ i * ncol + j];
			if( check == 0.0 )
				ok &= X_x[i] == 0.0;
			else
				ok &= std::fabs( X_x[i] / check - 1.0 ) <= eps;
		}

		// done with this X
		cholmod_free_dense(&X, &com);
	}

	// free memory
	cholmod_free_triplet(&T,    &com);
	cholmod_free_sparse( &H,    &com);
	cholmod_free_factor( &L,    &com);
	cholmod_free_dense(  &B,    &com);

	// finish up
	cholmod_finish(&com);

	// check count of malloc'ed - free'd
	ok &= com.malloc_count == 0;

	return ok;
}
// END C++