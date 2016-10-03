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
$begin sparse_ad_cholesky_ctor$$
$spell
	cholesky cholesky
	CppAD
	Alow
	const
	Eigen
$$

$section Sparse AD Cholesky Constructor$$

$head Syntax$$
$codei%CppAD::mixed::sparse_ad_cholesky cholesky(%Alow%)%$$

$head Public / Private$$
This is a public member function of the class $code sparse_ad_cholesky$$.
On the other hand, this class,
and all of its members, are implementation details and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Alow$$
This matrix has prototype
$codei%
	const Eigen::SparseMatrix<double, Eigen::ColMajor>& %Alow%
%$$
and is the lower triangle of a positive definite matrix.
Only positive definite matrices with the same sparsity pattern
can be factored using the $icode cholesky$$ object; i.e.,
square matrices with the same column size and
same set of possibly non-zero entries.

$end
*/
sparse_ad_cholesky::sparse_ad_cholesky(const sparse_d_matrix& Alow)
: CppAD::atomic_base<double> (
	"sparse_ad_cholesky",
	CppAD::atomic_base<double>::set_sparsity_enum
),
nc_( size_t( Alow.cols() ) )
{	assert( Alow.rows() == Alow.cols() );
	//
	// initialize ok flag as true
	// 2DO: change this to an error message similar to ipopt_fixed class
	ok_ = true;
	//
	// Step 1: Set Alow_pattern_
	assert( Alow_pattern_.row.size() == 0 );
	CppAD::mixed::sparse_eigen2info(Alow, Alow_pattern_);
	// Alow_pattern_.row and Alow_pattern_.col do not change
	// Alow_pattern_.val.size() does not change
	//
	// Step 2: analyze the sparsity pattern
	ldlt_obj_.analyzePattern( Alow );
	// This is the only call to ldlt_obj_.analyzePattern
	//
	// Step 3: Compute the Cholesky factor for this Alow
	// and the permutation all Alow (used with this object).
	ldlt_obj_.factorize( Alow );
	ok_ &= ldlt_obj_.info() == Eigen::Success;
	//
	// Step 4: Retrieve the permutation P_;
	P_   = ldlt_obj_.permutationP();
	//
	// Step 5: Set L_pattern_
	sparse_d_matrix L  =  ldlt_obj_.matrixL();
	assert( L_pattern_.row.size() == 0 );
	CppAD::mixed::sparse_eigen2info(L, L_pattern_);
	// L_pattern_.row and L_pattern_.col do not change
	// L_pattern_.val.size() does not change
}
/*
------------------------------------------------------------------------------
$begin sparse_ad_cholesy_ad$$
$spell
	Cholesky
	Alow
	CppAD
	const
	Eigen
$$

$section Using Sparse AD Cholesky Operation$$

$head Syntax$$
$icode%cholesky%.ad(%aAlow%, %aL%)%$$

$head Public / Private$$
This is a public member function of the class $code sparse_ad_cholesky$$.
On the other hand, this class,
and all of its members, are implementation details and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head aALow$$
This matrix has prototype
$codei%
	const Eigen::SparseMatrix< CppAD::AD<double>, Eigen::ColMajor>& %aAlow%
%$$
and is the lower triangle of a positive definite matrix.
It must have the same sparsity pattern as
$cref/Alow/sparse_ad_cholesky_ctor/Alow/$$ in the $icode cholesky$$
constructor.

$head aL$$
This matrix has prototype
$codei%
	Eigen::SparseMatrix< CppAD::AD<double>, Eigen::ColMajor>& %aL%
%$$
The input value of its elements does not matter.
Upon return, it is a lower triangular matrix such that
$latex \[
	P \cdot A \cdot P^\R{T} = L \cdot L^\R{T}
\]$$
where $latex L$$ is the lower triangular matrix corresponding to
$icode aL$$,
$latex P$$ is the permutation matrix corresponding to
$icode%
	%P% = %cholesky%.permutation()
%$$
and $latex A$$ is the symmetric positive definite matrix
with lower triangle equal to $icode aAlow$$.

$children%example/private/sparse_ad_cholesky_ad.cpp
%$$
$head Example$$
The file $cref sparse_ad_cholesky_ad.cpp$$ is an example
and test using this operation.

$end
*/
void sparse_ad_cholesky::ad(
	const sparse_ad_matrix& aAlow  ,
	sparse_ad_matrix&       aL     )
{	assert( nc_ == size_t( aAlow.rows() ) );
	assert( nc_ == size_t( aAlow.cols() ) );
	// -----------------------------------------------------------
	// packed version of Alow
	size_t nx = Alow_pattern_.row.size();
	CppAD::vector< CppAD::AD<double> > ax( nx );
	size_t index = 0;
	for(size_t j = 0; j < nc_; j++)
	{	for(sparse_ad_matrix::InnerIterator itr(aAlow, j); itr; ++itr)
		{	assert( Alow_pattern_.row[index] == size_t( itr.row() ) );
			assert( Alow_pattern_.col[index] == size_t( itr.col() ) );
			ax[ index ] = itr.value();
			++index;
		}
	}
	assert( index == nx );
	// -------------------------------------------------------------------
	// make call to packed vector verison of the atomic function
	size_t ny = L_pattern_.row.size();
	CppAD::vector< CppAD::AD<double> > ay( ny );
	(*this)(ax, ay);
	// -------------------------------------------------------------------
	// unpack ay into aL
	aL.resize(nc_, nc_);
	index = 0;
	for(size_t ell = 0; ell < ny; ell++)
	{	size_t i = L_pattern_.row[ell];
		size_t j = L_pattern_.col[ell];
		aL.insert(i, j) = ay[ index++ ];
	}
	assert( index == ny );
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
	{	for(size_t ell = 0; ell < nx; ell++)
			Alow_pattern_.val[ell] = tx[ ell * n_order + k ];
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
	// This is a very dumb algorithm that over estimates which elements
	// of L are variables. 2DO: create a much better estimate
	bool var = false;
	for(size_t j = 0; j < nx; j++)
		var |= vx[j];
	for(size_t i = 0; i < ny; i++)
		vy[i] = var;
	//
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
	{	for(size_t ell = 0; ell < nx; ell++)
			Alow_pattern_.val[ell] = tx[ ell * n_order + k ];
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
		for(size_t ell = 0; ell < nx; ell++)
			px[ ell * n_order + k ] = Alow_pattern_.val[ell];
	}
	//
	return true;
}

} } // END_CPPAD_MIXED_NAMESPASE

