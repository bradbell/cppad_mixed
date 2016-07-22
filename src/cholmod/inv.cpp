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
$begin ldlt_cholmod_inv$$

$section Compute a Subset of the Inverse of Factored Matrix$$
$spell
	ldlt_obj
	inv
	CppAD
	const
	cholmod
$$

$head Under Construction$$

$head Syntax$$
$icode%ldlt_obj%.inv(%row%, %col%, %val%)%$$

$head Prototype$$
$srcfile%src/cholmod/inv.cpp
	%0%// BEGIN_PROTOTYPE%// END_PROTOTYPE%1%$$

$head Private$$
The $cref ldlt_cholmod$$ class is an
$cref/implementation detail/ldlt_cholmod/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
This function solves for a subset of the inverse of the
sparse symmetric matrix that has been factored.

$head ldlt_obj$$
This object has prototype
$codei%
	const CppAD::mixed::ldlt_cholmod %ldlt_obj%
%$$
In addition, it must have a previous call to
$cref ldlt_cholmod_update$$.

$head row$$
This vector contains the row indices for the components of the
inverse that we are computing.

$head col$$
This vector contains the column indices for the components of the
inverse that we are computing. It must have the same size as $icode row$$.

$head val$$
This matrix must have the same size as $icode row$$.
The input values of its elements do not matter.
Upon return, it
contains the values for the components of the inverse that we are computing.
To be specific, for $icode%k% = 0 , %...%, %K%-1%$$,
$icode%val%[%k%]%$$
is $icode%row%[%k%]%$$, $icode%col%[%k%]%$$ component of the inverse.

$end
------------------------------------------------------------------------------
*/
# include <cstddef>
# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cppad/mixed/sparseinv.hpp>
# include <cppad/utility/index_sort.hpp>

namespace CppAD { namespace mixed {

// BEGIN_PROTOTYPE
void ldlt_cholmod::inv(
	const CppAD::vector<size_t>& row    ,
	const CppAD::vector<size_t>& col    ,
	CppAD::vector<double>&       val    )
// END_PROTOTYPE
{	using CppAD::vector;
	//
	// 2DO: move some of this work to the init routine.
	//
	assert( col.size() == row.size() );
	assert( val.size() == row.size() );
	//
	assert( factor_->n        == nrow_ );
	assert( factor_->minor    == nrow_ );
	assert( factor_->is_ll    == CHOLMOD_FALSE );
	assert( factor_->is_super == CHOLMOD_FALSE );
	assert( factor_->xtype    == CHOLMOD_REAL );
	//
	// information in the factor
	int*    L_perm  = reinterpret_cast<int*>( factor_->Perm );
	int*    L_p     = reinterpret_cast<int*>( factor_->p );
	int*    L_i     = reinterpret_cast<int*>( factor_->i );
	double* L_x     = reinterpret_cast<double*>( factor_->x );
	//
	// check order assumptions on factor
# ifndef NDEBUG
	for(int j = 0; j < int( nrow_) ; j++)
	{	// diagonal entry first
		assert( L_i[ L_p[j] ] == j );
		for(int k = L_p[j] + 1; k < L_p[j+1]; k++)
		{	// for each column, rows are sorted
			assert( L_i[k-1] < L_i[k] );
		}
	}
# endif
	//
	// determine column major order for matrix with rows and columns permuted
	size_t K = row.size();
	vector<size_t> key(K), ind(K);
	for(size_t k = 0; k < K; k++)
	{	size_t i = static_cast<size_t>( L_perm[ row[k] ] );
		size_t j = static_cast<size_t>( L_perm[ col[k] ] );
		key[k]   = i + j * nrow_;
	}
	CppAD::index_sort(key, ind);
	//
	// number of entries in use in factor
	size_t nnz  = static_cast<size_t>( L_p[nrow_] );
	//
	// number of rows and columns in the factor
	ptrdiff_t n  = static_cast<ptrdiff_t>(nrow_);
	//
	// convert L_p from int to ptrdiff_t
	vector<ptrdiff_t> Lp(nrow_ + 1);
	for(size_t j = 0; j <= nrow_; j++)
		Lp[j] = static_cast<ptrdiff_t>( L_p[j] );
	//
	// convert L_i from int to ptrdiff_t
	vector<ptrdiff_t> Li(nnz);
	for(size_t k = 0; k < nnz; k++)
		Li[k] = static_cast<ptrdiff_t>( L_i[k] );
	//
	// remove diagonal from L_x
	vector<double> Lx(nnz);
	for(size_t k = 0; k < nnz; k++)
		Lx[k] = L_x[k];
	//
	// extract the diagonal of the factor
	vector<double> d(nrow_);
	for(size_t j = 0; j < nrow_; j++)
	{	d[j] = Lx[ Lp[j] ];
		assert( d[j] != 0.0 );
		Lx[ Lp[j] ] = 0.0;
	}
	//
	// convert row, col to cholmod sparsity pattern
	vector<ptrdiff_t> Zp(nrow_+1), Zi(K);
	size_t k = 0;
	Zp[0]    = k;
	for(size_t j = 0; j < nrow_; j++)
	{	bool more = k < K;
		while( more )
		{	size_t c = static_cast<size_t>( L_perm[ col[ ind[k] ] ] );
			assert( c >= j );
			more = c == j;
			if( more )
			{	Zi[k] = static_cast<size_t>( L_perm[ row[ ind[k] ] ] );
				k++;
			}
			more = k < K;
		}
		Zp[j+1] = k;
	}
	// place where inverse is returned
	vector<double> Zx(K);
	//
	// work space
	vector<double> z(nrow_);
	vector<ptrdiff_t> Zdiagp(nrow_), Lmunch(nrow_);
	//
	// compute the subset of the inverse of the permuted matrix
	sparseinv(
		n,
		Lp.data(),
		Li.data(),
		Lx.data(),
		d.data(),
		Lp.data(),
		Li.data(),
		Lx.data(),
		Zp.data(),
		Zi.data(),
		Zx.data(),
		z.data(),
		Zdiagp.data(),
		Lmunch.data()
	);
	//
	// compute the inverse permutation
	vector<size_t> inv_perm(nrow_);
	for(size_t k = 0; k < nrow_; k++)
		inv_perm[ L_perm[k] ] = k;
	//
	// return the values
	for(size_t k = 0; k < K; k++)
		val[ ind[k] ] = Zx[ inv_perm[k] ];
	//
	return;
}

} }
