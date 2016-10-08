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
$begin sparse_ad_cholesy_ad$$
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

$head aALow$$
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
	size_t cj  = 0; // L index in column major order
	size_t rj  = 0; // L index in row major order
	for(size_t j = 0; j < nc_; j++)
	{	// Determine which elements of column j of L are variables using
		// B_{i,j} = sum_k L_{i,k} * L_{j,k}
		// where B = P * A * P^T and note that
		//
		// advance rj to the beginning of j-th row of L
		while( L_pattern_.row[ L_row_major_[rj] ] < j )
			++rj;
		//
		// first element in column j of L must be L(j,j)
		assert( L_pattern_.row[cj] == j );
		assert( L_pattern_.col[cj] == j );
		size_t Ljj = cj;
		//
		// There must be an element in row j of L
		assert( L_pattern_.row[ L_row_major_[rj] ] == j );
		//
		size_t ri = rj; // initialize ri to the beginning of row j
		while( cj < ny && L_pattern_.col[cj] == j )
		{	// This row index in L
			size_t i = L_pattern_.row[cj];
			//
			// Advance to beginning of row in in L
			while(
				L_pattern_.col[ L_row_major_[ri] ] <= j &&
				L_pattern_.row[ L_row_major_[ri] ] < i  )
				++ri;
			//
			// Determine of L(i,j) is a variable using
			// B(i,j) = L(i,0)*L(j,0) + ... + L(i,j)*L(j,j).
			// where B = P * A * P^T
			bool var = false;
			//
			// Check if element B(i,j) is a variable
			size_t r;
			size_t c;
			bool flag = true;
			while( flag )
			{	r = p_indices[ Alow_pattern_.row[ Alow_permuted_[ib] ] ];
				c = p_indices[ Alow_pattern_.col[ Alow_permuted_[ib] ] ];
				if( c > r )
					std::swap(r, c);
				flag = c < j || (c == j && r < i);
				if( flag )
					++ib;
			}
			if( c == j && r == i && vx[ Alow_permuted_[ib] ] )
			{	// Element B(i,j) is a variable
				var = true;
			}
			//
			size_t ck = 0; // index for column major access to L
			for(size_t k = 0; k <= j; k++)
			{	// check if L(i,k) * L(j,k) is a variable
				//
				// L(nc, nc) is last element in both row and column major
				// order so now worry about going off the end of the array
				while(
					L_pattern_.row[ L_row_major_[ri] ] <= i &&
					L_pattern_.col[ L_row_major_[ri] ] < k  )
					++ri;
				// found L(i,k)
				bool found_ik =
					L_pattern_.row[ L_row_major_[ri] ] == i &&
					L_pattern_.col[ L_row_major_[ri] ] == k;
				//
				while( L_pattern_.col[ck] <= k && L_pattern_.row[ck] <  i )
					++ck;
				// found L(j,k)
				bool found_jk =
					L_pattern_.col[ck] == k && L_pattern_.row[ck] == i;
				//
				if( found_ik && found_jk )
				{	if( k == j )
					{	if( i > j )
						{	assert( cj > Ljj );
							var |= vy[Ljj];
						}
					}
					else
					{	assert( L_row_major_[ri] < cj );
						assert( ck < cj );
						var |= vy[ L_row_major_[ri] ] || vy[ck];
					}
				}
			}
			vy[cj++] = var;
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

} } // END_CPPAD_MIXED_NAMESPASE

