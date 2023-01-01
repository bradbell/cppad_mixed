/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin sparseinv_xam.cpp$$
$spell
	sparseinv
	Cholmod
$$

$nospell
$bold This is cppad_mixed--20220519 documentation:$$ Here is a link to its
$href%https://cppad-mixed.readthedocs.io%current documentation%$$.
$$
$section Example Using sparseinv With Cholmod Factors$$

$head Example Description$$
This example uses a Cholmod LDLT factor to
compute the inversion of a
sparse symetric matrix $latex A$$
on the sets of indices where $latex A$$ is non-zero.
$srcthisfile%0%// BEGIN C++%// END C++%1%$$
$end
*/
// BEGIN C++
# include <cppad/mixed/include_cholmod.hpp>
# include <cppad/mixed/sparseinv.hpp>
# include <cppad/utility/sparse_rc.hpp>
# include <cppad/utility/vector.hpp>
# include <limits>
# include <cmath>
# include <cassert>
# include <vector>

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
	void matrix_prod(
		size_t n,
		bool A_transpose,
		const double* A,
		bool B_transpose,
		const double* B,
		double *C
	)
	{	for(size_t i = 0; i < n; ++i)
		{	for(size_t j = 0; j < n; ++j)
			{	C[i * n + j] = 0;
				for(size_t k = 0; k < n; ++k)
				{	double aik = A[i * n + k];
					if( A_transpose )
						aik = A[k * n + i];
					//
					double bkj = B[k * n + j];
					if( B_transpose )
						bkj = B[j * n + k];
					//
					C[i * n + j] += aik * bkj;
				}
			}
		}
		return;
	}
	/*
	---------------------------------------------------------------------------
	Debugging print routines
	void print_int(std::string name, size_t m, size_t n, const int* value)
	{	assert( m > 0 && n > 0 );
		using std::printf;
		//
		std::string indent = std::string( name.size() + 3, ' ');
		for(size_t i = 0; i < m; i++)
		{	if( i == 0 )
				printf("%s = ", name.c_str());
			else
				printf("%s", indent.c_str());
			for(size_t j = 0; j < n; ++j)
			{	if( j == 0 )
					printf("[");
				else
					printf(", ");
				printf("%d", value[i * n + j]);
			}
			printf("]\n");
		}
		return;
	}
	void print_double(std::string name, size_t m, size_t n, const double* value)
	{	assert( m > 0 && n > 0 );
		using std::printf;
		//
		std::string indent = std::string( name.size() + 3, ' ');
		for(size_t i = 0; i < m; i++)
		{	if( i == 0 )
				printf("%s = ", name.c_str());
			else
				printf("%s", indent.c_str());
			for(size_t j = 0; j < n; ++j)
			{	if( j == 0 )
					printf("[");
				else
					printf(", ");
				printf("%f", value[i * n + j]);
			}
			printf("]\n");
		}
		return;
	}
	---------------------------------------------------------------------------
	*/
}

bool sparseinv_xam(void)
{	bool ok      = true;
	double eps99 = 99.0 * std::numeric_limits<double>::epsilon();
	//
	// The matrix A
	size_t n = 4;
	const double A[] = {
		1.0, 0.0, 0.0,  1.0,
		0.0, 1.0, 0.0,  0.0,
		0.0, 0.0, 1.0,  0.0,
		1.0, 0.0, 0.0, 11.0
	};
	// The inverse of A is
	const double A_inv[] = {
		1.1, 0.0, 0.0, -0.1,
		0.0, 1.0, 0.0,  0.0,
		0.0, 0.0, 1.0,  0.0,
		-0.1, 0.0, 0.0, 0.1
	};
	//
	// check that I = A * A_inv  is the identity
	std::vector<double> I(n*n);
	matrix_prod(n, false, A, false, A_inv, I.data());
	for(size_t i = 0; i < n; ++i)
	{	for(size_t j = 0; j < n; ++j)
		{	double I_ij = double( i == j );
			ok &= ( std::fabs(I_ij - I[i * n + j]) < eps99 );
		}
	}
	//
	// initialize common
	cholmod_common com;
	cholmod_start(&com);
	// always do simplicial factorization
	com.supernodal = CHOLMOD_SIMPLICIAL;
	// do LDLT factorization and leave in LDLT form
	com.final_ll = CHOLMOD_FALSE;
	//
	// allocate triplet
	size_t nrow = n;
	size_t ncol = n;
	size_t nzmax = nrow * ncol;
	int    T_stype = CHOLMOD_STYPE_LOWER_TRIANGLE;
	int    T_xtype = CHOLMOD_REAL;
	cholmod_triplet* triplet =
		cholmod_allocate_triplet(nrow, ncol, nzmax, T_stype, T_xtype, &com);
	ok &= triplet->nnz ==  0;
	//
	// triplet and A_nnz
	size_t A_nnz = 0;
	for(size_t i = 0; i < nrow; i++)
	{	for(size_t j = 0; j <= i; j++)
		if( A[i * n + j] != 0.0 )
		{	add_T_entry(triplet, int(i), int(j), A[i * n + j] );
			++A_nnz;
			if( i != j )
				++A_nnz;
		}
	}
	//
	// convert triplet to sparse representation of A
	cholmod_sparse* sparse = cholmod_triplet_to_sparse(triplet, 0, &com);
	//
	//  the matrix
	cholmod_factor* factor = cholmod_analyze(sparse, &com);
	int flag               = cholmod_factorize(sparse, factor, &com);
	//
	// check properties of factor
	ok &= ( flag          == CHOLMOD_TRUE );  // return flag OK
	ok &= ( factor->n     == n );             // number of rows and coluns
	ok &= ( factor->minor == n );             // successful factorization
	ok &= ( factor->is_ll == CHOLMOD_FALSE ); // factorization is LDLT
	ok &= ( com.status == CHOLMOD_OK  );      // no problem with factorization
	//
	// information that defines the factor
	int*    L_perm = reinterpret_cast<int*>( factor->Perm );
	int*    L_p    = reinterpret_cast<int*>( factor->p );
	int*    L_i    = reinterpret_cast<int*>( factor->i );
	double* L_x    = reinterpret_cast<double*>( factor->x );
	//
	// P
	std::vector<double> P(n * n, 0.0);
	for(size_t i = 0; i < n; ++i)
		P[i * n + L_perm[i]] = 1.0;
	//
	// D
	std::vector<double> D(n);
	for(size_t j = 0; j < n; ++j)
		D[j] = L_x[ L_p[j] ];
	// ----------------------------------------------------------------------
	// Z_rc : sparsity patter for L*D*L^T + P*A*P^T
	// ----------------------------------------------------------------------
	//
	// PAPT_rc where PAPT = P*A*P^T
	CppAD::sparse_rc< CppAD::vector<size_t> > PAPT_rc(n, n, A_nnz);
	{	size_t k = 0;
		for(size_t i = 0; i < n; ++i)
		{	for(size_t j = 0; j < n; ++j)
			{	if( A[i * n + j] != 0.0 )
					PAPT_rc.set(k++, L_perm[i], L_perm[j]);
			}
		}
		ok &= (k == A_nnz);
	}
	// PAPT_col_major
	CppAD::vector<size_t> PAPT_col_major = PAPT_rc.col_major();
	//
	// Z_nnz
	size_t Z_nnz = A_nnz;
	{	size_t j = 0;
		int    m = L_p[j];
		size_t i = L_i[m];
		for(size_t k = 0; k < A_nnz; ++k)
		{	size_t r = PAPT_rc.row()[ PAPT_col_major[k] ];
			size_t c = PAPT_rc.col()[ PAPT_col_major[k] ];
			//
			// check if we need to advance m
			while( j < c || ( (j == c) && i <= r ) )
			{	// add (i, j) and (j, i) are in L * D * L^T sparsity
				// and not in P * A * P^T sparsity
				if( i != r && j != c )
				{	++Z_nnz;
					// if (i, j) not in PAPT_rc neither is (j, i)
					if( i != j )
						++Z_nnz;
				}
				// advance i, j, m
				++m;
				if( m == L_p[j+1] )
					++j;
				if( m < L_p[n] )
					i = L_i[m];
				else
				{	assert( j == n );
					i = n;
				}
			}
		}
	}
	//
	// Z_rc
	CppAD::sparse_rc< CppAD::vector<size_t> > Z_rc(n, n, Z_nnz);
	{	size_t j   = 0;
		int    m   = L_p[j];
		size_t i   = L_i[m];
		size_t ell = 0;
		for(size_t k = 0; k < A_nnz; ++k)
		{	size_t r = PAPT_rc.row()[ PAPT_col_major[k] ];
			size_t c = PAPT_rc.col()[ PAPT_col_major[k] ];
			Z_rc.set(ell++, r, c);
			//
			// check if we need to advance m
			while( j < c || ( (j == c) && i <= r ) )
			{	// add (i, j) and (j, i) are in L * D * L^T sparsity
				// and not in P * A * P^T sparsity
				if( i != r && j != c )
				{	Z_rc.set(ell++, i, j);
					// if (i, j) not in PAPT_rc neither is (j, i)
					if( i != j )
						Z_rc.set(ell++, j, i);
				}
				// advance i, j, m
				++m;
				if( m == L_p[j+1] )
					++j;
				if( m < L_p[n] )
					i = L_i[m];
				else
				{	assert( j == n );
					i = n;
				}
			}
		}
	}
	// ----------------------------------------------------------------------
	//
	// Z_p, Z_i
	CppAD::vector<size_t> Z_col_major = Z_rc.col_major();
	std::vector<int> Z_p(n+1), Z_i(Z_nnz);
	{	size_t k = 0;
		for(size_t j = 0; j < n; ++j)
		{	size_t ell = Z_col_major[k];
			ok &= (Z_rc.col()[ell] == j);
			Z_p[j] = int(k);
			Z_i[k] = int( Z_rc.row()[ell] );
			++k;
			while( k < Z_nnz && Z_rc.col()[ Z_col_major[k] ] == j )
			{	Z_i[k] = int( Z_rc.row()[ Z_col_major[k] ] );
				++k;
			}
			assert( k < Z_nnz || j == n-1 );
		}
		assert( k == Z_nnz );
		Z_p[n] = int( Z_nnz );
	}
	//
	// Z_x
	std::vector<double> Z_x(Z_nnz);
	//
	// work space
	std::vector<double> z(n, 0.0);
	std::vector<int> Zdiagp(n), Lmunch(n);
	sparseinv(
		static_cast<int>(n),
		L_p,
		L_i,
		L_x,
		D.data(),
		L_p,
		L_i,
		L_x,
		Z_p.data(),
		Z_i.data(),
		Z_x.data(),
		z.data(),
		Zdiagp.data(),
		Lmunch.data()
	);
	//
	// Z
	std::vector<double> Z(n * n, 0.0);
	for(size_t j = 0; j < n; ++j)
	{	for(size_t m = size_t(Z_p[j]); m < size_t(Z_p[j+1]); ++m)
		{	size_t i     = size_t( Z_i[m] );
			Z[i * n + j] = Z_x[m];
		}
	}
	//
	// check that P^T * Z * P = A_inv
	std::vector<double> ZP(n * n), PTZP(n * n);
	matrix_prod(n, false, Z.data(), false, P.data(), ZP.data());
	matrix_prod(n, true, P.data(), false, ZP.data(), PTZP.data());
# ifndef NDEBUG
	for(size_t k = 0; k < n * n; ++k)
	{	double abs_err = std::fabs( PTZP[k] - A_inv[k] );
		assert( abs_err < eps99 );
	}
# endif
	//
	// free memory
	cholmod_free_triplet(&triplet,    &com);
	cholmod_free_sparse( &sparse,    &com);
	cholmod_free_factor( &factor,    &com);
	//
	// finish up
	cholmod_finish(&com);
	//
	// check count of malloc'ed - free'd
	ok &= com.malloc_count == 0;
	//
	return ok;
}
// END C++
