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
$begin ldlt_cholmod_inv$$

$section Compute a Subset of the Inverse of Factored Matrix$$
$spell
	suitesparse
	suitesparse
	ldlt_obj
	inv
	CppAD
	const
	cholmod
	sparseinv
	nrow
	suitesparse
$$

$head Under Construction$$

$head Syntax$$
$icode%ldlt_obj%.inv(%row_in%, %col_in%, %val_out%)%$$

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

$head row_in$$
This vector contains the row indices for the components of the
inverse that we are computing.
It is assumed that these indices are the same for every call to
$code inv$$.

$head col_in$$
This vector contains the column indices for the components of the
inverse that we are computing. It must have the same size as $icode row_in$$.
It is assumed that these indices are the same for every call to
$code inv$$.

$head val_out$$
This matrix must have the same size as $icode row_in$$.
The input values of its elements do not matter.
Upon return, it
contains the values for the components of the inverse that we are computing.
To be specific, for $icode%k% = 0 , %...%, %K%-1%$$,
$icode%val_out%[%k%]%$$
is $icode%row_in%[%k%]%$$, $icode%col_in%[%k%]%$$ component of the inverse.

$head Method$$
This routine uses the
$cref/suitesparse/install_unix/System Requirements/suitesparse/$$ routine
$codei%
	MATLAB_Tools/sparseinv/sparseinv.c
%$$
to solve for the requested components of the inverse.

$head Member variables$$
The first time this routine is called, the following member variable
as set to their permanent values:

$subhead sparseinv_p_$$
This member variable has the prototype
$codei%
	CppAD::vector<int> sparseinv_p_
%$$
Its initial size is zero.
After the first call to $code inv$$, it had size $code nrow_+1$$.
For $icode%j% = 0, %...%, nrow_%$$,
$codei%sparseinv_p_[%j%]%$$ is where the row indices for the $th j$$
column start.

$subhead sparseinv_i_$$
This member variable has the prototype
$codei%
	CppAD::vector<int> sparseinv_i_
%$$
Its initial size is zero.
After the first call to $code inv$$, it had size
$codei
	sparseinv_p_[nrow_]
%$$
For each column index $icode j$$, and for
$codei
	%k% = sparseinv_p_[%j%] , %...%, sparseinv_p[%j%+1]-1
%$$
$codei%sparseinv_i_[%k%]%$$ is the index of the next entry in the $th j$$
column of the sparse matrix.
Note that for each column, the corresponding row indices are in order; i.e.
sorted.
This sparse matrix includes all the indices in
$icode row_in$$, $icode col_in$$,
all the indices in the factorization $code factor_$$,
and the transpose of the indices in $code factor_$$.

$subhead out2sparseinv_order_$$
This member variable has the prototype
$codei%
	CppAD::vector<size_t> out2sparseinv_order_
%$$
Its initial size is zero.
After the first call to $code inv$$,
it has the same size as $icode val_out$$.
It is a mapping
from the indices in $icode val_out$$ to the corresponding index in
the $code sparseinv_i$$ structures.
Note that all the indices in
$icode row_in$$ and $icode col_in$$ also appear in the representation
$code sparseinv_p_$$, $code sparseinv_i_$$.

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
	const CppAD::vector<size_t>& row_in    ,
	const CppAD::vector<size_t>& col_in    ,
	CppAD::vector<double>&       val_out   )
// END_PROTOTYPE
{	using CppAD::vector;
	//
	assert( col_in.size() == row_in.size() );
	assert( val_out.size() == row_in.size() );
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
	// number of rows and columns in the factor
	int n  = static_cast<int>(nrow_);
	//
	// extract the diagonal of the factor
	vector<double> d(nrow_);
	for(size_t j = 0; j < nrow_; j++)
	{	d[j] = L_x[ L_p[j] ];
		assert( d[j] != 0.0 );
	}
	//
	// number of requested indices and indices in L
	size_t K_in    = row_in.size();
	size_t K_L     = static_cast<size_t>( L_p[nrow_] );
	size_t K_total = K_in + 2 * K_L;
	// --------------------------------------------------------------------
	if( out2sparseinv_order_.size() == 0 )
	{
		// permuted version of requested indices
		vector<size_t> row(K_total), col(K_total);
		for(size_t k = 0; k < K_in; k++)
		{	row[k] = static_cast<size_t>( L_perm[ row_in[k] ] );
			col[k] = static_cast<size_t>( L_perm[ col_in[k] ] );
		}
		//
		// indices that are in L and L'
		for(size_t j = 0; j < nrow_; j++)
		{	for(int k = L_p[j]; k < L_p[j+1]; k++)
			{	size_t i = static_cast<size_t>( L_i[k] );
				// L
				row[K_in + 2 * k] = i;
				col[K_in + 2 * k] = j;
				// L'
				row[K_in + 2 * k + 1] = j;
				col[K_in + 2 * k + 1] = i;
			}
		}
		//
		// column major order
		vector<size_t> key(K_total), ind(K_total);
		for(size_t k = 0; k < K_total; k++)
			key[k] = row[k] + col[k] * nrow_;
		CppAD::index_sort(key, ind);
		//
		// put sparsity pattern as requred by sparseinv
		sparseinv_p_.resize(nrow_+1);
		sparseinv_i_.resize(0);
		// index conversize from sparse_mat_info to sparseinv format
		out2sparseinv_order_.resize(K_in);
		size_t k = 0;
		sparseinv_p_[0]    = int( sparseinv_i_.size() );
		for(size_t j = 0; j < nrow_; j++)
		{	bool more = k < K_total;
			while( more )
			{	size_t r = row[ ind[k] ];
				size_t c = col[ ind[k] ];
				assert( c >= j );
				more = c == j;
				if( more )
				{	sparseinv_i_.push_back( int(r) );
					if( ind[k] < K_in )
						out2sparseinv_order_[ ind[k] ] = sparseinv_i_.size()-1;
					k++;
					bool duplicate = k < K_total;
					while( duplicate )
					{	duplicate = (row[ind[k]] == r) & (col[ind[k]] == c);
						if( duplicate )
						{	if( ind[k] < K_in )
								out2sparseinv_order_[ ind[k] ] =
									sparseinv_i_.size()-1;
							k++;
							duplicate = k < K_total;
						}
					}
					more = k < K_total ;
				}
			}
			sparseinv_p_[j+1] = int( sparseinv_i_.size() );
		}
	}
	// ------------------------------------------------------------------------
	// place where inverse is returned
	vector<double> Zx( sparseinv_i_.size() );
	//
	// work space
	vector<double> z(nrow_);
	vector<int> Zdiagp(nrow_), Lmunch(nrow_);
	//
	// compute the subset of the inverse of the permuted matrix
	sparseinv(
		n,
		L_p,
		L_i,
		L_x,
		d.data(),
		L_p,
		L_i,
		L_x,
		sparseinv_p_.data(),
		sparseinv_i_.data(),
		Zx.data(),
		z.data(),
		Zdiagp.data(),
		Lmunch.data()
	);
	//
	// return the values
	for(size_t k = 0; k < K_in; k++)
		val_out[k] = Zx[ out2sparseinv_order_[k] ];
	//
	return;
}

} }
