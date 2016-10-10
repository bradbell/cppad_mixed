// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/sparse_ad_cholesky.hpp>
//
# include <cppad/mixed/sparse_low_tri_sol.hpp>
# include <cppad/mixed/sparse_up_tri_sol.hpp>
# include <cppad/mixed/sparse_scale_diag.hpp>
# include <cppad/mixed/sparse_low2sym.hpp>
# include <cppad/mixed/sparse_mat2low.hpp>
# include <cppad/mixed/sparse_eigen2info.hpp>
# include <cppad/mixed/sparse_info2eigen.hpp>
# include <iostream>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
// ============================================================================
// Public member functions
// ============================================================================
/*
$begin sparse_ad_cholesky_initialize$$
$spell
	cholesky cholesky
	CppAD
	Alow
	const
	Eigen
	Cholesky
$$

$section Initialize Sparse AD Cholesky Factorization$$

$head Syntax$$
$icode%cholesky%(%ad_Alow%)%$$

$head Public / Private$$
This is a public member function of the class $code sparse_ad_cholesky$$.
On the other hand, this class,
and all of its members, are implementation details and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head cholesky$$
This is object has prototype
$codei%
	sparse_ad_cholesky %cholesky%
%$$
and was created with the default constructor.
The $code initialize$$ routine should be called once
for each $icode cholesky$$ object.

$head ad_Alow$$
This matrix has prototype
$codei%
	const Eigen::SparseMatrix< CppAD::AD<double>, Eigen::ColMajor>& %ad_Alow%
%$$
and is the lower triangle of a positive definite matrix.
Only positive definite matrices
$cref/A/sparse_ad_cholesky/Notation/A/$$ with the same sparsity pattern
can be factored using the $icode cholesky$$ object; i.e.,
square matrices with the same column size and
same set of possibly non-zero entries.

$head Restriction$$
The $code CppAD::AD<double>$$ tape cannot be recording when this
function is called and hence all such AD objects are parameters
(not variables).

$end
*/
void sparse_ad_cholesky::initialize(const sparse_ad_matrix& ad_Alow)
{	assert( ad_Alow.rows() == ad_Alow.cols() );
	// ---------------------------------------------------------------------
	// number of rows and columns in Alow and  L
	nc_ = ad_Alow.rows();
	// ----------------------------------------------------------------------
	// double version of Alow
	sparse_d_matrix Alow(nc_, nc_);
	for(size_t j = 0; j < nc_; j++)
	{	for(sparse_ad_matrix::InnerIterator itr(ad_Alow, j); itr; ++itr)
		{	assert( Parameter( itr.value() ) );
			Alow.insert( itr.row(), itr.col() ) = Value( itr.value() );
		}
	}
	// ----------------------------------------------------------------------
	// initialize ok_ as true
	// 2DO: change this to an error message similar to ipopt_fixed class
	ok_ = true;
	// ----------------------------------------------------------------------
	// Set Alow_pattern_
	assert( Alow_pattern_.row.size() == 0 );
	CppAD::mixed::sparse_eigen2info(Alow, Alow_pattern_);
	// Alow_pattern_.row and Alow_pattern_.col do not change
	// Alow_pattern_.val.size() does not change
	// ----------------------------------------------------------------------
	// analyze the Alow sparsity pattern using ldlt_obj_
	ldlt_obj_.analyzePattern( Alow );
	// This is the only call to ldlt_obj_.analyzePattern
	// ----------------------------------------------------------------------
	// Compute the Cholesky factor for this Alow
	// and the permutation all Alow that are used with this object.
	ldlt_obj_.factorize( Alow );
	ok_ &= ldlt_obj_.info() == Eigen::Success;
	// ----------------------------------------------------------------------
	// Retrieve the permutation P_;
	P_   = ldlt_obj_.permutationP();
	// ----------------------------------------------------------------------
	// Set L_pattern_
	sparse_d_matrix L  =  ldlt_obj_.matrixL();
	assert( L_pattern_.row.size() == 0 );
	CppAD::mixed::sparse_eigen2info(L, L_pattern_);
	// L_pattern_.row and L_pattern_.col do not change
	// L_pattern_.val.size() does not change
	// ----------------------------------------------------------------------
	// a vector of integers corresponding to the permutation
	Eigen::Matrix<size_t, Eigen::Dynamic, 1> p_indices(nc_);
	for(size_t i = 0; i < nc_; i++)
		p_indices[i] = i;
	p_indices = P_ * p_indices;
	// ----------------------------------------------------------------------
	// Indices that sort lower triangle of P * A * P^T in column major order
	size_t nx = Alow_pattern_.row.size();
	Alow_permuted_.resize(nx);
	CppAD::vector<size_t> keys( nx );
	for(size_t ia = 0; ia < nx; ia++)
	{	size_t i = p_indices[ Alow_pattern_.row[ia] ];
		size_t j = p_indices[ Alow_pattern_.col[ia] ];
		if( j > i )
			std::swap(i, j);
		keys[ia] = j * nc_ + i;
	}
	CppAD::index_sort(keys, Alow_permuted_);
	// ----------------------------------------------------------------------
	// Indices that sort L_pattern_ in row major order
	size_t ny = L_pattern_.row.size();
	L_row_major_.resize(ny);
	keys.resize( ny );
	for(size_t ell = 0; ell < ny; ell++)
	{	size_t i = L_pattern_.row[ell];
		size_t j = L_pattern_.col[ell];
		assert( i < nc_ );
		assert( j < nc_ );
		keys[ell] = i * nc_ + j;
	}
	CppAD::index_sort(keys, L_row_major_);
}
/*
$begin sparse_ad_cholesky_p$$
$spell
	Cholesky
	CppAD
$$

$section Using Sparse AD Cholesky Permutation P$$

$head Syntax$$
$icode%P% = %cholesky%.permutation()%$$

$head Prototype$$
$srcfile%src/eigen/sparse_ad_cholesky.cpp
	%4%// BEGIN PERMUTATION PROTOTYPE%// END PERMUTATION PROTOTYPE%1%$$

$head Public / Private$$
This is a public member function of the class $code sparse_ad_cholesky$$.
On the other hand, this class,
and all of its members, are implementation details and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head cholesky$$
This is object has prototype
$codei%
	sparse_ad_cholesky %cholesky%
%$$
The $cref/initialize/sparse_ad_cholesky_initialize/$$ routine must be called
before using $code permutation$$ for a $icode cholesky$$ object.

$head P$$
The return value is the permutation matrix
$cref/P/sparse_ad_cholesky/Notation/P/$$.
The permutation corresponding to $icode cholesky$$ does not change.

$head Example$$
The file $cref sparse_ad_chol_perm.cpp$$ is an example
and test using this operation.

$end
*/
// BEGIN PERMUTATION PROTOTYPE
const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>&
sparse_ad_cholesky::permutation(void)
// END PERMUTATION PROTOTYPE
{	return P_; }
/*
------------------------------------------------------------------------------
$begin sparse_ad_cholesky_eval$$
$spell
	Cholesky
	Alow
	CppAD
	const
	Eigen
	eval
$$

$section Using Sparse AD Cholesky Factor L$$

$head Syntax$$
$icode%cholesky%.eval(%ad_Alow%, %ad_L%)%$$

$head Public / Private$$
This is a public member function of the class $code sparse_ad_cholesky$$.
On the other hand, this class,
and all of its members, are implementation details and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head ad_Alow$$
This matrix has prototype
$codei%
	const Eigen::SparseMatrix< CppAD::AD<double>, Eigen::ColMajor>& %ad_Alow%
%$$
and is the lower triangle of a positive definite matrix
$cref/A/sparse_ad_cholesky/Notation/A/$$.
It must have the same sparsity pattern as
$cref/ad_Alow/sparse_ad_cholesky_initialize/ad_Alow/$$ in the $icode cholesky$$
initialization.

$head ad_L$$
This matrix has prototype
$codei%
	Eigen::SparseMatrix< CppAD::AD<double>, Eigen::ColMajor>& %ad_L%
%$$
The input value of its elements does not matter.
Upon return, it is a lower triangular matrix
$cref/L/sparse_ad_cholesky/Notation/L/$$ corresponding to the matrix $latex A$$
specified by $icode ad_Alow$$.

$head Example$$
The file $cref sparse_ad_chol_eval.cpp$$ is an example
and test using this operation.

$end
*/
void sparse_ad_cholesky::eval(
	const sparse_ad_matrix& ad_Alow  ,
	sparse_ad_matrix&       ad_L     )
{	assert( nc_ == size_t( ad_Alow.rows() ) );
	assert( nc_ == size_t( ad_Alow.cols() ) );
	// -------------------------------------------------------------------
	// packed version of Alow
	size_t nx = Alow_pattern_.row.size();
	CppAD::vector< CppAD::AD<double> > ax( nx );
	size_t ia = 0;
	for(size_t j = 0; j < nc_; j++)
	{	for(sparse_ad_matrix::InnerIterator itr(ad_Alow, j); itr; ++itr)
		{	assert( Alow_pattern_.row[ia] == size_t( itr.row() ) );
			assert( Alow_pattern_.col[ia] == size_t( itr.col() ) );
			ax[ ia ] = itr.value();
			++ia;
		}
	}
	assert( ia == nx );
	// -------------------------------------------------------------------
	// make call to packed vector verison of the atomic function
	size_t ny = L_pattern_.row.size();
	CppAD::vector< CppAD::AD<double> > ay( ny );
	(*this)(ax, ay);
	// -------------------------------------------------------------------
	// unpack ay into ad_L
	ad_L.resize(nc_, nc_);
	for(size_t ell = 0; ell < ny; ell++)
	{	size_t i = L_pattern_.row[ell];
		size_t j = L_pattern_.col[ell];
		ad_L.insert(i, j) = ay[ell];
	}
	return;
}
/*
==============================================================================
private functions
==============================================================================

$begin set_jac_sparsity$$
$spell
	Jacobian
	jac
	cholesky
	CppAD
	Alow
$$

$section Set the Jacobian Sparsity Pattern$$

$head Syntax$$
$codei%set_jac_sparsity(%jac_sparsity%)%$$

$head Private$$
This is a private member function of the class $code sparse_ad_cholesky$$.
In addition, this class,
and all of its members, are implementation details and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head jac_sparsity$$
This argument has one of the following prototypes
$codei%
	CppAD::sparse_list& %jac_sparsity%
	CppAD::sparse_pack& %jac_sparsity%
%$$
Its input value does not matter.
Upon return it contains the jacobian sparsity pattern for the mapping from
the lower triangle of the symmetric matrix
$cref/Alow/sparse_ad_cholesky_eval/ad_Alow/$$ to the Cholesky factor
$cref/L/sparse_ad_cholesky_eval/ad_L/$$.

$subhead n_set$$
The number of sets $icode%jac_sparsity%.n_set()%$$ is equal
to the number of elements in the Cholesky factor $latex L$$.

$subhead end$$
The end for each set $icode%jac_sparsity%.end()%$$ is equal
to the number of elements in the lower triangle of the symmetric
matrix $latex A$$. The elements of the set are greater than or equal zero
and less than the end value.

$subhead elements$$
If $icode j$$ is an element of set $icode i$$,
then the $th i$$ element of $latex L$$ depends on the
$th j$$ element of the lower triangle of $latex A$$
(using column major ordering).


$end
*/
template <class Sparsity>
void sparse_ad_cholesky::set_jac_sparsity(Sparsity& jac_sparsity)
{
	// -------------------------------------------------------------------
	// Initialize Jacobian sparsity pattern for L as empty
	size_t nx = Alow_pattern_.row.size();
	size_t ny = L_pattern_.row.size();
	jac_sparsity.resize(ny, nx);
	//
	// a vector of integers corresponding to the permutation
	Eigen::Matrix<size_t, Eigen::Dynamic, 1> p_indices(nc_);
	for(size_t i = 0; i < nc_; i++)
		p_indices[i] = i;
	p_indices = P_ * p_indices;
	//
	// Determine sparsity pattern for L
	size_t ib  = 0; // Blow index in column major order
	size_t cij = 0; // index of L(i,j) in column major order
	size_t rj  = 0; // index of row L(j,:) in row major order
	for(size_t j = 0; j < nc_; j++)
	{	// Determine sparsity pattern for j-th column of L
		//
		// advance rj to the beginning of j-th row of L
		while( L_pattern_.row[ L_row_major_[rj] ] < j )
			++rj;
		//
		// first element in column j of L must be L(j,j)
		assert( L_pattern_.row[cij] == j );
		assert( L_pattern_.col[cij] == j );
		//
		// There must be an element in row j of L
		assert( L_pattern_.row[ L_row_major_[rj] ] == j );
		//
		size_t ri = rj; // initialize ri to the beginning of row j
		while( cij < ny && L_pattern_.col[cij] == j )
		{	// The row index in L corresponding to cij
			size_t i = L_pattern_.row[cij];
			//
			// Advance ri to beginning of row i in L
			while(
				L_pattern_.col[ L_row_major_[ri] ] <= j &&
				L_pattern_.row[ L_row_major_[ri] ] < i  )
				++ri;
			//
			// Determine sparsity pattern for L(i,j) using
			// B(i,j) = L(i,0) * L(j,0) + ... + L(i,j) * L(j,j).
			// where B = P * A * P^T
			//
			// Check for element B(i,j)
			size_t r, c;
			bool flag = true;
			while( flag )
			{	// advance ib to next element
				r = p_indices[ Alow_pattern_.row[ Alow_permuted_[ib] ] ];
				c = p_indices[ Alow_pattern_.col[ Alow_permuted_[ib] ] ];
				if( c > r )
					std::swap(r, c);
				flag = c < j || (c == j && r < i);
				if( flag )
					++ib;
			}
			// check if next element in B is B(i,j)
			if( r == i && c == j )
			{	// L(i,j) depends on element of Alow that corresponds to B(i,j)
				jac_sparsity.add_element(cij, Alow_permuted_[ib] );
			}
			//
			size_t rik = ri; // initialize index for next element in row i
			size_t rjk = rj; // initialize index for next element in row j
			for(size_t k = 0; k <= j; k++)
			{	// check if L(i,k) * L(j,k) is a variable
				//
				// L(nc, nc) is last element in both row and column major
				// order so now worry about going off the end of the array
				assert( L_pattern_.row[ L_row_major_[rik] ] >= i );
				assert( L_pattern_.row[ L_row_major_[rjk] ] >= j );
				//
				while(
					L_pattern_.row[ L_row_major_[rik] ] == i &&
					L_pattern_.col[ L_row_major_[rik] ] < k  )
					++rik;
				// found L(i,k)
				bool found_ik =
					L_pattern_.row[ L_row_major_[rik] ] == i &&
					L_pattern_.col[ L_row_major_[rik] ] == k;
				//
				while(
					L_pattern_.row[ L_row_major_[rjk] ] == j &&
					L_pattern_.col[ L_row_major_[rjk] ] < k  )
					++rjk;
				// found L(j,k)
				bool found_jk =
					L_pattern_.row[ L_row_major_[rjk] ] == j &&
					L_pattern_.col[ L_row_major_[rjk] ] == k;
				//
				if( found_ik && found_jk )
				{	// both L(i,k) and L(j,k) can be non-zero
					size_t cik = L_row_major_[rik];
					size_t cjk = L_row_major_[rjk];
					if( k == j )
					{	if( i > j )
						{	// L(i,j) * L(j,j)
							assert( cjk < cij );
							// set cij = set cij union set cjk
							jac_sparsity.binary_union(
								cij,
								cij,
								cjk,
								jac_sparsity
							);
						}
						else
						{	// L(j,j) * L(j,j) case is determining cjk
							assert( cjk == cij );
						}
					}
					else
					{	// L(i,k) * L(j,k)
						assert( cik < cij );
						assert( cjk < cij );
						//
						// set cij = set cij union set cik
						jac_sparsity.binary_union(
							cij,
							cij,
							cik,
							jac_sparsity
						);
						// set cij = set cij union set cjk
						jac_sparsity.binary_union(
							cij,
							cij,
							cjk,
							jac_sparsity
						);
					}
				}
			}
			// Done with sparsity for L(i,j)
			cij++;
		}
	}
	return;
}
// ===========================================================================
// These virtual functions are specified by the CppAD::atomic_base requirements
// ===========================================================================
// forward mode routine for this operation
bool sparse_ad_cholesky::forward(
	// lowest order Taylor coefficient we are evaluating
	size_t                          p ,
	// highest order Taylor coefficient we are evaluating
	size_t                          q ,
	// which components of x are variables
	const CppAD::vector<bool>&      vx ,
	// which components of y are variables
	CppAD::vector<bool>&            vy ,
	// tx [ j * (q+1) + k ] is x_j^k
	const CppAD::vector<double>&    tx ,
	// ty [ i * (q+1) + k ] is y_i^k
	CppAD::vector<double>&          ty )
{	//
	assert( p <= q );
	//
	size_t n_order = q + 1;
	size_t nx      = Alow_pattern_.row.size();
	size_t ny      = L_pattern_.row.size();
	assert( vx.size() == 0 || nx == vx.size() );
	assert( vx.size() == 0 || ny == vy.size() );
	assert( nx * n_order == tx.size() );
	assert( ny * n_order == ty.size() );
	// -------------------------------------------------------------------
	// f_Alow and f_L
	CppAD::vector<sparse_d_matrix> f_Alow(n_order), f_L(n_order);
	// -------------------------------------------------------------------
	// unpack tx into f_Alow
	for(size_t k = 0; k < n_order; k++)
	{	for(size_t ia = 0; ia < nx; ia++)
			Alow_pattern_.val[ia] = tx[ ia * n_order + k ];
		CppAD::mixed::sparse_info2eigen(f_Alow[k], Alow_pattern_, nc_, nc_);
	}
	// -------------------------------------------------------------------
	// for orders less than p, unpack ty into f_L
	for(size_t k = 0; k < p; k++)
	{	// unpack f_L values for this order
		for(size_t ell = 0; ell < ny; ell++)
			L_pattern_.val[ell] = ty[ ell * n_order + k];
		CppAD::mixed::sparse_info2eigen(f_L[k], L_pattern_, nc_, nc_);
	}
	if( p == 0 )
	{	// compute result for zero order
		ldlt_obj_.factorize( f_Alow[0] );
		if( ldlt_obj_.info() != Eigen::Success )
		{	ok_ = false;
			return false;
		}
		Eigen::VectorXd D2 = sqrt( ldlt_obj_.vectorD().array() ).matrix();
		f_L[0]            =  ldlt_obj_.matrixL() * D2.asDiagonal();
		assert( P_.indices() == ldlt_obj_.permutationP().indices() );
	}
	Eigen::SparseMatrix<double, Eigen::RowMajor> L0 = f_L[0];
	// compute result for orders greater than or equal 1 and p.
	for(size_t k = std::max(p, size_t(1)); k < n_order; k++)
	{	// convert Alow to A_k
		sparse_d_matrix A_k = CppAD::mixed::sparse_low2sym(f_Alow[k]);
		//
		// initialize sum as E_k =  P * A_k * P^T
		sparse_d_matrix tmp1  = P_ * A_k;
		sparse_d_matrix f_sum = tmp1 * P_.transpose();
		//
		// compute E_k - B_k
		for(size_t ell = 1; ell < k; ell++)
			f_sum -= f_L[ell] * f_L[k-ell].transpose();
		//
		// compute L_0^{-1} * (E_k - B_k)
		tmp1 = CppAD::mixed::sparse_low_tri_sol(L0, f_sum);
		//
		// compute L_0^{-1} * (E_k - B_k) * L_0^{-T}
		sparse_d_matrix tmp2 = CppAD::mixed::sparse_low_tri_sol(
			L0, tmp1.transpose()
		).transpose();
		//
		// divide the diagonal by 2
		CppAD::mixed::sparse_scale_diag(0.5, tmp2);
		//
		// low[ L_0^{-1} * (E_k - B_k) * L_0^{-T} ]
		// L_k = L_0 * low[  L_0^{-1} * (E_k - B_k) * L_0^{-T} ]
		f_L[k] = L0 * CppAD::mixed::sparse_mat2low(tmp2);
	}
	// -------------------------------------------------------------------
	// pack f_L into ty
	// -------------------------------------------------------------------
	for(size_t k = p; k < n_order; k++)
	{	CppAD::mixed::sparse_eigen2info(f_L[k], L_pattern_);
		for(size_t ell = 0; ell < ny; ell++)
			ty[ ell * n_order + k ] = L_pattern_.val[ell];
	}
	// -------------------------------------------------------------------
	// check if we are computing vy
	if( vx.size() == 0 )
		return true;
	// -------------------------------------------------------------------
	// a vector of integers corresponding to the permutation
	Eigen::Matrix<size_t, Eigen::Dynamic, 1> p_indices(nc_);
	for(size_t i = 0; i < nc_; i++)
		p_indices[i] = i;
	p_indices = P_ * p_indices;
	// -------------------------------------------------------------------
	// Determine which components of L are variables.
	//
	// Determine which elemements of L are variables in column major order
	size_t ib  = 0; // Blow index in column major order
	size_t cij = 0; // index of L(i,j) in column major order
	size_t rj  = 0; // index of row L(j,:) in row major order
	for(size_t j = 0; j < nc_; j++)
	{	// Determine which elements of j-th column of L are variables
		//
		// advance rj to the beginning of j-th row of L
		while( L_pattern_.row[ L_row_major_[rj] ] < j )
			++rj;
		//
		// first element in column j of L must be L(j,j)
		assert( L_pattern_.row[cij] == j );
		assert( L_pattern_.col[cij] == j );
		//
		// There must be an element in row j of L
		assert( L_pattern_.row[ L_row_major_[rj] ] == j );
		//
		size_t ri = rj; // initialize ri to the beginning of row j
		while( cij < ny && L_pattern_.col[cij] == j )
		{	// The row index in L corresponding to cij
			size_t i = L_pattern_.row[cij];
			//
			// Advance ri to beginning of row i in L
			while(
				L_pattern_.col[ L_row_major_[ri] ] <= j &&
				L_pattern_.row[ L_row_major_[ri] ] < i  )
				++ri;
			//
			// Determine of L(i,j) is a variable using
			// B(i,j) = L(i,0) * L(j,0) + ... + L(i,j) * L(j,j).
			// where B = P * A * P^T
			bool var = false;
			//
			// Check if element B(i,j) is a variable
			size_t r, c;
			bool flag = true;
			while( flag )
			{	// advance ib to next element
				r = p_indices[ Alow_pattern_.row[ Alow_permuted_[ib] ] ];
				c = p_indices[ Alow_pattern_.col[ Alow_permuted_[ib] ] ];
				if( c > r )
					std::swap(r, c);
				flag = c < j || (c == j && r < i);
				if( flag )
					++ib;
			}
			// check if next element in B is B(i,j) and it is a variable
			if( c == j && r == i && vx[ Alow_permuted_[ib] ] )
			{	// Element L(i,j) depends on B(i,j) which is a variable
				var = true;
			}
			//
			size_t rik = ri; // initialize index for next element in row i
			size_t rjk = rj; // initialize index for next element in row j
			for(size_t k = 0; k <= j; k++)
			{	// check if L(i,k) * L(j,k) is a variable
				//
				// L(nc, nc) is last element in both row and column major
				// order so now worry about going off the end of the array
				assert( L_pattern_.row[ L_row_major_[rik] ] >= i );
				assert( L_pattern_.row[ L_row_major_[rjk] ] >= j );
				//
				while(
					L_pattern_.row[ L_row_major_[rik] ] == i &&
					L_pattern_.col[ L_row_major_[rik] ] < k  )
					++rik;
				// found L(i,k)
				bool found_ik =
					L_pattern_.row[ L_row_major_[rik] ] == i &&
					L_pattern_.col[ L_row_major_[rik] ] == k;
				//
				while(
					L_pattern_.row[ L_row_major_[rjk] ] == j &&
					L_pattern_.col[ L_row_major_[rjk] ] < k  )
					++rjk;
				// found L(j,k)
				bool found_jk =
					L_pattern_.row[ L_row_major_[rjk] ] == j &&
					L_pattern_.col[ L_row_major_[rjk] ] == k;
				//
				if( found_ik && found_jk )
				{	// both L(i,k) and L(j,k) can be non-zero
					size_t cik = L_row_major_[rik];
					size_t cjk = L_row_major_[rjk];
					if( k == j )
					{	if( i > j )
						{	// L(i,j) * L(j,j)
							// Element L(i,j) depends on L(j,j)
							assert( cjk < cij );
							var |= vy[cjk];
						}
						else
						{	// L(j,j) * L(j,j) case is determining cjk
							assert( cjk == cij );
						}
					}
					else
					{	// L(i,k) * L(j,k)
						assert( cik < cij );
						assert( cjk < cij );
						//
						// L(i,j) depends on L(i,k)
						var |= vy[ cik ];
						// L(i,j) depends on L(j,k)
						var |= vy[ cjk ];
					}
				}
			}
			// set variable value for L(i,j)
			vy[cij++] = var;
		}
	}
	return true;
}
// ------------------------------------------------------------------
// reverse mode routine for this operation
bool sparse_ad_cholesky::reverse(
	// highest order Taylor coefficient that we are computing derivative of
	size_t                     q ,
	// forward mode Taylor coefficients for x variables
	const CppAD::vector<double>&     tx ,
	// forward mode Taylor coefficients for y variables
	const CppAD::vector<double>&     ty ,
	// upon return, derivative of G[ F[ {x_j^k} ] ] w.r.t {x_j^k}
	CppAD::vector<double>&           px ,
	// derivative of G[ {y_i^k} ] w.r.t. {y_i^k}
	const CppAD::vector<double>&     py )
{	//
	size_t n_order = q + 1;
	size_t nx      = Alow_pattern_.row.size();
	size_t ny      = L_pattern_.row.size();
	//
	assert( nx * n_order == tx.size() );
	assert( ny * n_order == ty.size() );
	assert( px.size()    == tx.size() );
	assert( py.size()    == ty.size() );
	//
	// -------------------------------------------------------------------
	// declare f_Alow, f_L, r_Alow, r_L
	CppAD::vector<sparse_d_matrix> f_Alow(n_order), f_L(n_order);
	CppAD::vector<sparse_d_matrix> r_Alow(n_order), r_L(n_order);
	// -------------------------------------------------------------------
	// unpack tx into f_Alow
	for(size_t k = 0; k < n_order; k++)
	{	for(size_t ia = 0; ia < nx; ia++)
			Alow_pattern_.val[ia] = tx[ ia * n_order + k ];
		CppAD::mixed::sparse_info2eigen(f_Alow[k], Alow_pattern_, nc_, nc_);
	}
	// -------------------------------------------------------------------
	// unpack ty into f_L
	for(size_t k = 0; k < n_order; k++)
	{	// unpack f_L values for this order
		for(size_t ell = 0; ell < ny; ell++)
			L_pattern_.val[ell] = ty[ ell * n_order + k];
		CppAD::mixed::sparse_info2eigen(f_L[k], L_pattern_, nc_, nc_);
	}
	// -------------------------------------------------------------------
	// unpack py into r_L
	for(size_t k = 0; k < n_order; k++)
	{	for(size_t ell = 0; ell < ny; ell++)
			L_pattern_.val[ell] = py[ ell * n_order + k];
		CppAD::mixed::sparse_info2eigen(r_L[k], L_pattern_, nc_, nc_);
	}
	// -------------------------------------------------------------------
	// initialize r_Alow as zero
	for(size_t k = 0; k < n_order; k++)
		r_Alow[k].resize(nc_, nc_);
	// -------------------------------------------------------------------
	// Cholesky factorization
	Eigen::SparseMatrix<double, Eigen::RowMajor> L0 = f_L[0];
	//
	// start at highest order and go down
	for(size_t k1 = n_order; k1 > 1; k1--)
	{	size_t k = k1 - 1;
		// L_0^T * bar{L}_k
		sparse_d_matrix tmp1 = L0.transpose() * r_L[k];
		//
		// low[ L_0^T * bar{L}_k ]
		CppAD::mixed::sparse_scale_diag(0.5, tmp1);
		sparse_d_matrix tmp2 = CppAD::mixed::sparse_mat2low(tmp1);
		//
		// L_0^{-T} * low[ L_0^T * bar{L}_k ]
		tmp1 = CppAD::mixed::sparse_up_tri_sol(L0.transpose(), tmp2);
		//
		// L_0^{-T} * low[ L_0^T * bar{L}_k ]^T * L_0^{-1}
		sparse_d_matrix Mk = CppAD::mixed::sparse_up_tri_sol(
			L0.transpose(), tmp1.transpose()
		).transpose();
		//
		// remove Lk, \bar{Alow}_k += P^T * M0 * P
		sparse_d_matrix barB_k = CppAD::mixed::sparse_low2sym(Mk);
		tmp1                   = P_.transpose() * barB_k;
		r_Alow[k]             += tmp1 * P_;
		//
		// compute barB_k
		barB_k   = - barB_k;
		//
		// remove C_k using
		// 2 * lower[ bar{B}_k L_k ]
		tmp1    = barB_k * f_L[k];
		r_L[0] += 2.0 * CppAD::mixed::sparse_mat2low( tmp1 );
		//
		// remove B_k
		for(size_t ell = 1; ell < k; ell++)
		{	// bar{L}_ell = 2 * lower( bar{B}_k * L_{k-ell} )
			tmp1      = barB_k * f_L[k-ell];
			r_L[ell] += 2.0 * CppAD::mixed::sparse_mat2low( tmp1 );
		}
	}
	// L_0^T * bar{L}_0
	sparse_d_matrix tmp1 = L0.transpose() * r_L[0];
	//
	// low[ L_0^T * bar{L}_0 ]
	CppAD::mixed::sparse_scale_diag(0.5, tmp1);
	sparse_d_matrix tmp2 = CppAD::mixed::sparse_mat2low( tmp1 );
	//
	// L_0^{-T} low[ L_0^T * bar{L}_0 ]
	tmp1 = CppAD::mixed::sparse_up_tri_sol(L0.transpose(), tmp2);
	//
	// M_0 = L_0^{-T} low[ L_0^T * bar{L}_0 ]^T L_0^{-1}
	sparse_d_matrix M0 = CppAD::mixed::sparse_up_tri_sol(
		L0.transpose(), tmp1.transpose()
	);
	// remove L0, \bar{Alow}_0 += 2.0 * low[ P^T * M0 * P ]
	tmp1 = P_.transpose() * M0;
	tmp2 = tmp1 * P_;
	CppAD::mixed::sparse_scale_diag(0.5, tmp2);
	r_Alow[0] += 2.0 * CppAD::mixed::sparse_mat2low( tmp2 );
	// ------------------------------------------------------------------
	// pack r_Alow into px
	for(size_t k = 0; k < n_order; k++)
	{	CppAD::mixed::sparse_eigen2info(r_Alow[k], Alow_pattern_);
		for(size_t ia = 0; ia < nx; ia++)
			px[ ia * n_order + k ] = Alow_pattern_.val[ia];
	}
	//
	return true;
}
// ---------------------------------------------------------------------------
// vectorBool forward Jacobian sparsity pattern for S = f'(x) * R
bool sparse_ad_cholesky::for_sparse_jac(
	// number of columns in the matrix R
	size_t                        q         ,
	// sparsity pattern for R
	const CppAD::vectorBool&      r         ,
	// sparsity pattern for S
	CppAD::vectorBool&            s         ,
	// parameters in argument to atomic function
	const CppAD::vector<double>&  not_used  )
{	// make sure we have boolean version of sparsity for f'(x)
	if( jac_sparsity_bool_.n_set() == 0 )
		set_jac_sparsity(jac_sparsity_bool_);
	//
	// number of elements in domain and range of f(x)
	size_t nx = Alow_pattern_.row.size();
	size_t ny = L_pattern_.row.size();
	//
	assert( r.size() == nx * q );
	assert( s.size() == ny * q );
	//
	// compute sparsity pattern for S = f'(x) * R
	for(size_t i = 0; i < ny; i++)
	{	for(size_t j = 0; j < q; j++)
		{	// initialize sparsity pattern for S(i,j)
			bool s_ij = false;
			//
			// S(i, j) = sum_k J(i, k) R(k, j) where J = f'(x)
			for(size_t k = 0; k < nx; k++)
			{	bool J_ik = jac_sparsity_bool_.is_element(i, k);
				bool R_kj = r[ k * q + j ];
				s_ij     |= (J_ik & R_kj);
			}
			// set sparsity pattern for S(i, j)
			s[ i * q + j ] = s_ij;
		}
	}
	return true;
}
// ---------------------------------------------------------------------------
// vectorBool reverse Jacobian sparsity pattern for S = R * f'(x)
bool sparse_ad_cholesky::rev_sparse_jac(
	// number of rows in the matrix R
	size_t                        q         ,
	// sparsity pattern for R^T
	const CppAD::vectorBool&      rt        ,
	// sparsity pattern for S^T
	CppAD::vectorBool&            st        ,
	// parameters in argument to atomic function
	const CppAD::vector<double>&  not_used  )
{	// make sure we have boolean version of sparsity for f'(x)
	if( jac_sparsity_bool_.n_set() == 0 )
		set_jac_sparsity(jac_sparsity_bool_);
	//
	// number of elements in domain and range of f(x)
	size_t nx = Alow_pattern_.row.size();
	size_t ny = L_pattern_.row.size();
	//
	assert( rt.size() == ny * q );
	assert( st.size() == nx * q );
	//
	// compute sparsity pattern for S^T = f'(x)^T * R^T
	for(size_t i = 0; i < q; i++)
	{	for(size_t j = 0; j < nx; j++)
		{	// initialize sparsity pattern for S(i,j)
			bool s_ij = false;
			//
			// S(i, j) = sum_k R(i, k) J(k, j) where J = f'(x)
			for(size_t k = 0; k < ny; k++)
			{	bool R_ik = rt[ k * q + i ];
				bool J_kj = jac_sparsity_bool_.is_element(k, j);
				s_ij     |= (R_ik & J_kj);
			}
			// set sparsity pattern for S^T(j, i)
			st[ j * q + i ] = s_ij;
		}
	}
	return true;
}
// ----------------------------------------------------------------------------
// vectorBool reverse Hessian sparsity for V(x) = (g o f)^(2) (x) * R
bool sparse_ad_cholesky::rev_sparse_hes(
	// variable flag for x arguments
	const vector<bool>&                   vx       ,
	// sparsity pattern for scalar valued S(x) = g'[ f(x) ]
	const vector<bool>&                   s        ,
	// sparsity pattern for T(x) = (g o f)' (x) = S(x) * f'(x)
	vector<bool>&                         t        ,
	// number of columns in the matrix R
	size_t                                q        ,
	// sparsity pattern for the matrix R
	const CppAD::vectorBool&              r        ,
	// sparsity pattern for U(x) = g^(2)[ f(x) ] f'(x) R
	const CppAD::vectorBool&              u        ,
	// sparsity pattern for V(x) = (g o f)^(2) (x) * R
	CppAD::vectorBool&                    v        ,
	// parameters in argument to atomic function
	const vector<double>&                 not_used )
{	// make sure we have boolean version of sparsity for f'(x)
	if( jac_sparsity_bool_.n_set() == 0 )
		set_jac_sparsity(jac_sparsity_bool_);
	//
	// number of elements in domain and range of f(x)
	size_t nx = Alow_pattern_.row.size();
	size_t ny = L_pattern_.row.size();
	//
	assert( vx.size() == nx );
	assert( s.size()  == ny );
	assert( t.size()  == nx );
	assert( r.size()  == nx * q );
	assert( u.size()  == ny * q );
	assert( v.size()  == nx * q );
	assert( jac_sparsity_bool_.n_set() == ny );
	assert( jac_sparsity_bool_.end()   == nx );
	//
	// compute atomic sparsity pattern for T(x) = S(x) * f'(x)
	for(size_t j = 0; j < nx; j++)
	{	// T(j) = sum_k S(i) * J(i, j)
		bool t_j = false;
		for(size_t k = 0; k < ny; k++)
			t_j |= ( s[k] & jac_sparsity_bool_.is_element(k, j) );
		t[j] = t_j;
	}
	//
	// not yet finished
	return false;
}

} } // END_CPPAD_MIXED_NAMESPASE

