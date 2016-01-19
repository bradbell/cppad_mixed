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
$begin cholmod_xam.cpp$$
$spell
	Cholmod
	Cholesky
$$

$section Example Using Cholmod Cholesky Factorization$$

$head Problem Description$$
We are given the matrix
$latex \[
	A = \left( \begin{array}{ccc}
		5 & 4 & 2 \\
		4 & 5 & 1 \\
		2 & 1 & 5
	\end{array} \right)
\] $$
We use $latex A^k$$ to denote the upper-left $latex k \times k$$
principal minor. The determinant of its principal minors are:
$latex \[
\begin{array}{rcl}
	\det \left( A^1 \right) & = & 5 \\
	\det \left( A^2 \right) & = & 9 \\
	\det \left( A^3 \right) & = & 36
\end{array}
\] $$
In addition
$latex \[
	A^{-1} = \frac{1}{36}
	\left( \begin{array}{ccc}
		24  & -18 & -6 \\
		-18 & 21  & 3 \\
		-6  & 3   & 9
	\end{array} \right)
\] $$
which can be checked by multiplying by $latex A$$.

$head Source Code$$
$code
$verbatim%example/cholmod_xam.cpp%5%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cholmod.h>
# include <limits>
# include <cmath>
# include <cassert>

# define CHOLMOD_TRUE  1
# define CHOLMOD_FALSE 0

namespace {
	void add_T_entry(cholmod_triplet *T, int r, int c, double x)
	{	size_t k = T->nnz;

		((int*)T->i)[k] = r;
		((int*)T->j)[k] = c;
		((double*)T->x)[k] = x;

		T->nnz++;
	}
}

bool cholmod_xam(void)
{	bool ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();

	double A_inv[] = {
		 24.0, -18.0, -6.0,
		-18.0,  21.0,  3.0,
		 -6.0,   3.0,  9.0
	};
	for(size_t i = 0; i < sizeof(A_inv)/sizeof(A_inv[0]); i++)
		A_inv[i] /= 36.;


	cholmod_common com;
	cholmod_start(&com);

	// always do simplicial factorization
	com.supernodal = CHOLMOD_SIMPLICIAL;

	// do LDL' factorization and leave in LDL' form
	com.final_ll = CHOLMOD_FALSE;


	// allocate triplet
	size_t nrow = 3;
	size_t ncol = 3;
	size_t nzmax = nrow * ncol;
	int    stype = -1;    // lower triangle specified for symmetric case
	int    xtype = CHOLMOD_REAL;
	cholmod_triplet *T =
		cholmod_allocate_triplet(nrow, ncol, nzmax, stype, xtype, &com);
	ok &= T->nnz ==  0;

	// triplet entries corresponding to lower triangle of A
	add_T_entry(T, 0, 0, 5.);
	add_T_entry(T, 1, 0, 4.);
	add_T_entry(T, 2, 0, 2.);
	add_T_entry(T, 1, 1, 5.);
	add_T_entry(T, 2, 1, 1.);
	add_T_entry(T, 2, 2, 5.);

	// convert triplet to sparse representation of A
	cholmod_sparse* A =
		cholmod_triplet_to_sparse(T, 0, &com);

	// factor the matrix
	cholmod_factor *L = cholmod_analyze(A, &com);
	cholmod_factorize(A, L, &com);

	// check properties of factor
	assert( L->n     == nrow );          // number of rows and coluns
	assert( L->minor == nrow );          // successful factorization
	assert( L->is_ll == CHOLMOD_FALSE ); // factorization is LDL'

	// compute log of determinant of diagonal D
	double log_det_A;
	int*    L_p  = (int *) L->p;
	int*    L_i  = (int *) L->i;
	int*    L_nz = (int *) L->nz;
	double* L_x  = (double *) L->x;
	for(size_t j = 0; j < nrow; j++)
	{	assert( L_nz[j] > 0 );
		bool found = false;
		for(int k = 0; k < L_nz[j]; k++)
		{	if( size_t( L_i[ L_p[j] + k ] ) == j )
			{	assert( ! found );
				found = true;
				// j-th element on diagonal of factorization
				double dj = L_x[ L_p[j] + k ];
				assert( dj > 0.0 );
				log_det_A += std::log(dj);
			}
			assert(found);
		}
	}
	// check its value
	ok &= std::fabs( log_det_A / std::log(36.0) - 1.0 ) <= eps;


	// sparsity pattern for right hand side column vector
	size_t Bset_nrow       = nrow;
	size_t Bset_ncol       = 1;
	size_t Bset_nzmax      = 1;
	int    Bset_sorted     = CHOLMOD_TRUE;
	int    Bset_packed     = CHOLMOD_TRUE;
	int    Bset_stype      = 0;
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

	// sparsity pattern for solution column vector
	size_t Xset_nzmax = nrow;
	cholmod_sparse* Xset = cholmod_allocate_sparse(
		Bset_nrow,
		Bset_ncol,
		Xset_nzmax,
		Bset_sorted,
		Bset_packed,
		Bset_stype,
		Bset_xtype,
		&com
	);

	// non-zero values in the identity matrix are always one.
	cholmod_dense *B = cholmod_ones(nrow, 1, xtype, &com);

	// place where the solution is returned
	cholmod_dense *X = cholmod_zeros(nrow, 1, xtype, &com);
	double *X_x = (double *) X->x;

	// work space vectors that can be reused
	cholmod_dense *Y = NULL;
	cholmod_dense *E = NULL;

	// Recover the lower triangle of the symmetrix inverse matrix
	for(size_t j = 0; j < ncol; j++)
	{	// j-th column of identity matrix
		int* Bset_p = (int *) Bset->p;
		int* Bset_i = (int *) Bset->i;
		Bset_p[0] = 0;   // column index
		Bset_p[1] = 1;   // number of non-zeros in column vector
		Bset_i[0] = j;   // row index

		// entire column vector
		int* Xset_p = (int *) Xset->p;
		int* Xset_i = (int *) Xset->i;
		Xset_p[0] = 0;        // column index
		Xset_p[1] = nrow - j; // number of non-zeros in column vector
		for(size_t i = j; i < nrow; i++)
			Xset_i[i-j] = (int) i; // row index

		// solve an (i, j) element of the inverse of A
		cholmod_solve2(
			CHOLMOD_A,
			L,
			B,
			Bset,
			&X,
			&Xset,
			&Y,
			&E,
			&com
		);
		size_t X_len = (size_t) Xset_p[1];
		ok &= X_len == nrow - j;
		for(size_t i = j; i < nrow; i++)
		{	ok &= Xset_i[i-j] == (int) i;
			ok &= std::fabs( X_x[i] / A_inv[i*ncol+j] - 1.0 ) <= eps;
		}
	}

	// free memory
	cholmod_free_sparse(&Xset, &com);
	cholmod_free_sparse(&Bset, &com);
	cholmod_free_dense(&E, &com);
	cholmod_free_dense(&Y, &com);
	cholmod_free_dense(&X, &com);
	cholmod_free_dense(&B, &com);
	cholmod_free_factor(&L, &com);
	cholmod_free_sparse(&A, &com);
	cholmod_free_triplet(&T, &com);

	// finish up
	cholmod_finish(&com);

	// check memory count
	ok &= com.memory_inuse == 0;

	return ok;
}
// END C++
