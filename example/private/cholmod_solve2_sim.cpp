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
$begin cholmod_solve2_sim.cpp$$
$spell
	Cholmod
	Cholesky
	Xset
	Bset
	sys
	Pt
$$

$section Cholmod Posterior Simulations Using Sparse Hessian of Likelihood$$

$head Purpose$$
Suppose that $latex A$$ is a sparse positive definite Hessian of
a likelihood at the maximum likelihood estimate for its unknown parameters.
The corresponding asymptotic covariance for posterior distribution
of the parameters is normal with covariance $latex A^{1}$$.
A vector $latex x$$ with this covariance can be simulated as
$latex \[
	x = R b
\] $$
where $latex R$$ is defined by $latex A^{-1} = R R^\R{T}$$ and
$latex b$$ is a normal with mean zero and the identity covariance.
Suppose we have used $code cholmod$$ to compute the factorization
$latex \[
	L D L^\R{T} = P A P^\R{T}
\] $$
where $latex L$$ is lower triangular, $latex D$$ is diagonal,
and $latex P$$ is a permutation matrix.
It follows that
$latex \[
	A = P^\R{T} L D L^\R{T} P
\] $$
$latex \[
	A^{-1} = P^\R{T} L^{-\R{T}} D^{-1} L^{-1} P
\] $$
$latex \[
	R = P^\R{T} L^{-\R{T}} D^{-1/2}
\] $$
$latex \[
	x = P^\R{T} L^{-\R{T}} D^{-1/2} b
\] $$
A vector $latex x$$ with covariance $latex A^{-1}$$
can be computed by using this formula where $latex b$$
is simulated with the identity for its covariance.

$head Example$$
Solve for this example
$latex \[
	A = \left( \begin{array}{ccc}
		1 & 1 & 0 \\
		1 & 5 & 0 \\
		0 & 0 & 4
	\end{array} \right)
\] $$
The corresponding cholesky factorization is
$latex \[
	L = \left( \begin{array}{ccc}
		1 & 0 & 0 \\
		2 & 1 & 0 \\
		0 & 0 & 1
	\end{array} \right)
	\W{,}
	D = \left( \begin{array}{ccc}
		1 & 0 & 0 \\
		0 & 1 & 0 \\
		0 & 0 & 2
	\end{array} \right)
	\W{,}
	P = \left( \begin{array}{ccc}
		1 & 0 & 0 \\
		0 & 1 & 0 \\
		0 & 0 & 1
	\end{array} \right)
\] $$
$latex \[
	L^{-1} = \left( \begin{array}{ccc}
		1  & 0 & 0 \\
		-2 & 1 & 0 \\
		0  & 0 & 1
	\end{array} \right)
	\W{,}
	L^{-T} = \left( \begin{array}{ccc}
		1 & -2 & 0 \\
		0 & 1  & 0 \\
		0 & 0  & 1
	\end{array} \right)
\] $$
It follows that, for this example
$latex \[
	x = L^{-\R{T}} ( b_0 , b_1 , b_2 / 2 )^\R{T}
\] $$
$latex \[
	x = ( b_0 - 2 b_1 , b_1 , b_2 / 2 )^\R{T}
\] $$

$head Source Code$$
$code
$srcfile%example/private/cholmod_solve2_sim.cpp%5%// BEGIN C++%// END C++%1%$$
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

bool cholmod_solve2_sim(void)
{	bool ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();

	// initialize common
	cholmod_common com;
	cholmod_start(&com);

	// always do simplicial factorization
	com.supernodal = CHOLMOD_SIMPLICIAL;

	// do LDL' factorization and leave in LDL' form
	com.final_ll = CHOLMOD_FALSE;

	// allocate triplet
	size_t ncol = 3;
	size_t nrow = ncol;
	size_t nzmax = nrow * ncol;
	int    T_stype = CHOLMOD_STYPE_LOWER_TRIANGLE;
	int    T_xtype = CHOLMOD_REAL;
	cholmod_triplet *T =
		cholmod_allocate_triplet(nrow, ncol, nzmax, T_stype, T_xtype, &com);
	ok &= T->nnz ==  0;

	// triplet entries corresponding to lower triangle of A
	add_T_entry(T, 0, 0, 1.); // A_00 = 1
	add_T_entry(T, 1, 0, 2.); // A_10 = 2
	add_T_entry(T, 1, 1, 5.); // A_11 = 5
	add_T_entry(T, 2, 2, 4.); // A_22 = 4
	ok &= T->nnz ==  4;

	// convert triplet to sparse representation of A
	cholmod_sparse* A = cholmod_triplet_to_sparse(T, 0, &com);

	// factor the matrix
	cholmod_factor *L = cholmod_analyze(A, &com);
# ifndef NDEBUG
	int flag =
# endif
	cholmod_factorize(A, L, &com);

	// pointer the the in index vectors for the factor
	int* L_p    = (int *)    L->p;
	double* L_x = (double *) L->x;
# ifndef NDEBUG
	int* L_i    = (int *)    L->i;
# endif

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

	// We will solve for all components of the right hand side
	int* Bset_p = (int *) Bset->p;
	int* Bset_i = (int *) Bset->i;
	Bset_p[0]   = 0;
	Bset_p[1]   = nrow;
	for(size_t i = 0; i < nrow; i++)
		Bset_i[i] = i;

	// set the vector B = b
	cholmod_dense *B = cholmod_zeros(nrow, 1, T_xtype, &com);
	double* B_x     = (double *) B->x;
	for(size_t i = 0; i < nrow; i++)
		B_x[i] = double( i * i + 1);

	// check for solution
	double check[3];
	check[0] = B_x[0] - 2.0 * B_x[1];
	check[1] = B_x[1];
	check[2] = B_x[2] / 2.0;

	// work space vectors that can be reused
	cholmod_dense *Y = NULL;
	cholmod_dense *E = NULL;

	// place where the solution is returned
	cholmod_dense *X = NULL;
	cholmod_sparse* Xset = NULL;

	// set B = D^{-1/2} b
	for(size_t j = 0; j < ncol; j++)
	{	// first element for each column is alwasy the diagonal element
		assert( size_t( L_i [ L_p[j] ] ) == j );
		// j-th element on the diagonal of D
		double dj = L_x[ L_p[j] ];
		assert( dj > 0.0 );
		B_x[j] = B_x[j] / std::sqrt( dj );
	}

	// set X = L^-T * D^{-1/2} b
	int sys = CHOLMOD_Lt;
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
	assert( ni == nrow );
	//
	// set B = L^-T * D^{-1/2} b
	for(size_t k = 0; k < ni; k++)
	{	size_t i = Xset_i[k];
		B_x[i]   = X_x[i];
	}

	// set X = P^T * L^-T * D^{-1/2} b
	sys = CHOLMOD_Pt;
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
	// The four arrays Xset, X, Y, E may change each time
	X_x     = (double *) X->x;
	Xset_p  = (int *) Xset->p;
	Xset_i  = (int *) Xset->i;
	ni      = size_t( Xset_p[1] );
	assert( ni == nrow );
	//
	for(size_t k = 0; k < ni; k++)
	{	size_t i = Xset_i[k];
		double xi =  X_x[i];
		ok &= std::fabs( xi / check[i] - 1.0 ) <= eps;
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
