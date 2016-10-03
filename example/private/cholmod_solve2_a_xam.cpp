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
$begin cholmod_solve2_a_xam.cpp$$
$spell
	Cholmod
	Cholesky
	Xset
	Bset
	sys
	Pt
$$

$section Example Using cholmod_solve2 Cholesky Factorization$$

$head Problem Description$$
Solve for $latex x$$ in the equation $latex A x = b$$ where $latex A$$
is defined below and $latex b$$ is a column of the identity matrix.
Hence the solution $latex x$$ is the corresponding column of
$latex A^{-1}$$.
$latex \[
	G = \left( \begin{array}{ccc}
		5 & 4 & 2 \\
		4 & 5 & 1 \\
		2 & 1 & 5
	\end{array} \right)
	\W{,}
	A = \left( \begin{array}{cc}
		G & 0 \\
		0 & G
	\end{array} \right)
\] $$
We use $latex G^k$$ to denote the upper-left $latex k \times k$$
principal minor. The determinant of its principal minors are:
$latex \[
	\det \left( G^1 \right) = 5  \W{,}
	\det \left( G^2 \right) = 9  \W{,}
	\det \left( G^3 \right) = 36
\] $$
Hence, $latex G$$ and $latex A$$ are positive definite.
In addition
$latex \[
	G^{-1} = \frac{1}{36}
	\left( \begin{array}{ccc}
		24  & -18 & -6 \\
		-18 & 21  & 3 \\
		-6  & 3   & 9
	\end{array} \right)
	\W{,}
	A^{-1} = \left( \begin{array}{cc}
		G^{-1} & 0 \\
		0 & G^{-1}
	\end{array} \right)
\] $$
which can be checked by multiplying by $latex G G^{-1}$$.

$head Bset$$
I think the $code cholmod_solve2$$ documentation would be clearer if
the sentence
'The entries in $code Bset$$ are a subset of $code Xset$$
(except if sys is $code CHOLMOD_P$$ or $code CHOLMOD_Pt$$ ).'
were changed to
'The entries in $code Bset$$ are a subset of $code Xset$$
and only these entries are assured to appear in $code Xset$$
(except if sys is $code CHOLMOD_P$$ or $code CHOLMOD_Pt$$ ).'

$head Source Code$$
$code
$srcfile%example/private/cholmod_solve2_a_xam.cpp%5%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cholmod.h>
# include <limits>
# include <cmath>
# include <cassert>

# define CHOLMOD_TRUE                  1
# define CHOLMOD_FALSE                 0
# define CHOLMOD_STYPE_NOT_SYMMETRIC   0
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

bool cholmod_solve2_a_xam(void)
{	bool ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();

	// Example of checking for cholmod-3.0.4
	// (which is distributed as part of SuiteSparse-4.4.3)
	// int version[3];
	// cholmod_version(version);
	// ok &= version[0] == 3;
	// ok &= version[1] == 0;
	// ok &= version[2] == 4;

	double A_inv[] = {
		 24.0, -18.0, -6.0,   0.0,   0.0,  0.0,
		-18.0,  21.0,  3.0,   0.0,   0.0,  0.0,
		 -6.0,   3.0,  9.0,   0.0,   0.0,  0.0,
		  0.0,   0.0,  0.0,  24.0, -18.0, -6.0,
		  0.0,   0.0,  0.0, -18.0,  21.0,  3.0,
		  0.0,   0.0,  0.0,  -6.0,   3.0,  9.0
	};
	size_t nrow = 6;
	for(size_t i = 0; i < nrow * nrow; i++)
		A_inv[i] /= 36.;

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
	// first G block
	add_T_entry(T, 0, 0, 5.); // A_00 = 5
	add_T_entry(T, 1, 0, 4.); // A_10 = 4
	add_T_entry(T, 2, 0, 2.); // A_20 = 2
	add_T_entry(T, 1, 1, 5.); // A_11 = 5
	add_T_entry(T, 2, 1, 1.); // A_21 = 1
	add_T_entry(T, 2, 2, 5.); // A_22 = 5
	// second G block
	add_T_entry(T, 3, 3, 5.); // A_33 = 5
	add_T_entry(T, 4, 3, 4.); // A_43 = 4
	add_T_entry(T, 5, 3, 2.); // A_53 = 2
	add_T_entry(T, 4, 4, 5.); // A_44 = 5
	add_T_entry(T, 5, 4, 1.); // A_54 = 1
	add_T_entry(T, 5, 5, 5.); // A_55 = 5

	// convert triplet to sparse representation of A
	cholmod_sparse* A = cholmod_triplet_to_sparse(T, 0, &com);

	// factor the matrix
	cholmod_factor *L = cholmod_analyze(A, &com);
# ifndef NDEBUG
	int flag =
# endif
	cholmod_factorize(A, L, &com);

	// check properties of factor
	assert( flag     == CHOLMOD_TRUE );  // return flag OK
	assert( L->n     == nrow );          // number of rows and coluns
	assert( L->minor == nrow );          // successful factorization
	assert( L->is_ll == CHOLMOD_FALSE ); // factorization is LDL'
	assert( com.status == CHOLMOD_OK );  // no problem with factorization

	// sparsity pattern for right hand side column vector
	size_t Bset_nrow       = nrow;
	size_t Bset_ncol       = 1;
	size_t Bset_nzmax      = nrow;
	int    Bset_sorted     = CHOLMOD_FALSE;
	int    Bset_packed     = CHOLMOD_TRUE;
	int    Bset_stype      = CHOLMOD_STYPE_NOT_SYMMETRIC;
	int    Bset_xtype      = CHOLMOD_PATTERN;
	cholmod_sparse* Bset = cholmod_allocate_sparse(
		Bset_nrow,
		Bset_ncol,
		Bset_nzmax,
		Bset_sorted,
		Bset_packed,
		Bset_stype,
		Bset_xtype,
		&com
	);
	// We need to solve for three components of the right hand side
	int* Bset_p = (int *) Bset->p;
	int* Bset_i = (int *) Bset->i;
	Bset_p[0]   = 0;
	Bset_p[1]   = 3;

	// set B to a vector of ones
	cholmod_dense *B = cholmod_ones(nrow, 1, T_xtype, &com);
	double* B_x     = (double *) B->x;

	// work space vectors that can be reused
	cholmod_dense *Y = NULL;
	cholmod_dense *E = NULL;

	// place where the solution is returned
	cholmod_dense *X = NULL;
	cholmod_sparse* Xset = NULL;

	// one by one recover columns of A_inv
	for(size_t j = 0; j < ncol; j++)
	{
		// set B_x to j-th column of identity matrix
		for(size_t i = 0; i < 6; i++)
			B_x[i] = 0.0;
		B_x[j] = 1.0;

		if( j < 3 )
		{	// we need to solve for first three components of X
			Bset_i[0] = 0;
			Bset_i[1] = 1;
			Bset_i[2] = 2;

			// value of B_x in rows 3, 4, 5 do not matter
			for(size_t i = 3; i < 6; i++)
				B_x[i] = 2.0;
		}
		else
		{	// we need to solve for last three components of X
			Bset_i[0] = 3;
			Bset_i[1] = 4;
			Bset_i[2] = 5;

			// value of B_x in rows 0, 1, 2 do not matter
			for(size_t i = 0; i < 3; i++)
				B_x[i] = 2.0;
		}

		// solve A * x = b
		int sys = CHOLMOD_A;
# ifndef NDEBUG
		flag =
# endif
		cholmod_solve2(
			sys,
			L,
			B,
			Bset,
			&X,
			&Xset,
			&Y,
			&E,
			&com
		);
		//
		assert( flag == CHOLMOD_TRUE ); // return flag OK
		assert( Xset->nrow   == nrow );
		assert( Xset->ncol   == 1    );
		assert( Xset->xtype  == CHOLMOD_PATTERN );
		assert( Xset->packed == CHOLMOD_TRUE);
		//
		// The four arrays Xset, X, Y, E may change each time
		double* X_x     = (double *) X->x;
		int*    Xset_p  = (int *) Xset->p;
		int*    Xset_i  = (int *) Xset->i;
		size_t  ni      = size_t( Xset_p[1] );
		//
		// solution vector for this right hand side
		double x[nrow];
		for(size_t i = 0; i < nrow; i++)
			x[i] = 0.0;
		for(size_t k = 0; k < ni; k++)
		{	size_t i = Xset_i[k];
			x[i]     = X_x[i];
		}
		for(size_t i = 0; i < nrow; i++)
		{	double check_i = A_inv[ i * nrow + j];
			if( check_i == 0.0 )
				ok &= x[i] == 0.0;
			else
				ok &= std::fabs( x[i] / check_i - 1.0 ) <= eps;
		}
	}

	// free memory
	cholmod_free_triplet(&T,    &com);
	cholmod_free_sparse( &A,    &com);
	cholmod_free_factor( &L,    &com);
	cholmod_free_sparse( &Bset, &com);
	cholmod_free_sparse( &Xset, &com);
	cholmod_free_dense(  &B,    &com);
	cholmod_free_dense(  &Y,    &com);
	cholmod_free_dense(  &E,    &com);
	cholmod_free_dense(  &X,    &com);

	// finish up
	cholmod_finish(&com);

	// check count of malloc'ed - free'd
	ok &= com.malloc_count == 0;

	return ok;
}
// END C++
