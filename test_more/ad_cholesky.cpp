// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <Eigen/SparseCore>
# include <Eigen/SparseCholesky>
# include <cppad/cppad.hpp>
# include <cppad/mixed/sparse_mat_info.hpp>
# include <iostream>

namespace { // BEGIN_EMPTY_NAMESPACE
// -----------------------------------------------------------------
// define AD, Dynamic, some vector and matrix types
// -----------------------------------------------------------------
using CppAD::AD;
using Eigen::Dynamic;
typedef CppAD::vector<size_t>                             s_vector;
typedef CppAD::vector<double>                             d_vector;
typedef CppAD::vector< AD<double> >                       ad_vector;
//
typedef Eigen::PermutationMatrix<Dynamic, Dynamic>        perm_matrix;
typedef Eigen::Matrix<double, Dynamic, Dynamic>           dense_d_matrix;
typedef Eigen::SparseMatrix<double>                       sparse_d_matrix;
typedef Eigen::SparseMatrix< AD<double> >                 sparse_ad_matrix;
// -----------------------------------------------------------------
// convert sparse_ad_matrix -> sparse_d_matrix
// -----------------------------------------------------------------
sparse_d_matrix admat2mat(
	const sparse_ad_matrix& amat   ,
	const s_vector&         nnz   )
{	size_t nr = size_t( amat.rows() );
	size_t nc = size_t( amat.cols() );
	assert( size_t( nnz.size() ) == nc );
	//
	sparse_d_matrix mat(nr, nc);
	mat.reserve(nnz);
	for(size_t j = 0; j < nc; j++)
	{	for(sparse_ad_matrix::InnerIterator itr(amat, j); itr; ++itr)
		{	AD<double> av = itr.value();
			mat.insert( itr.row(), itr.col() ) = Value( av );
		}
	}
	return mat;
}
// -----------------------------------------------------------------
class atomic_ad_cholesky : public CppAD::atomic_base<double> {
public:
	// -----------------------------------------------------------------
	// data (in order that it is initialized)
	// -----------------------------------------------------------------
	// OK flag
	bool ok_;
	//
	// Number of non-zeros and sparsity pattern for Alow (set by constructor)
	s_vector                      Alow_nnz_;
	CppAD::mixed::sparse_mat_info Alow_pattern_;

	// Object used for Cholesky factorization
	// (analyzePattern by constructor only)
	Eigen::SimplicialLDLT<sparse_d_matrix> ldlt_obj_;

	// Value of the permutation matrix (set by constructor)
	perm_matrix P_;

	// Value of the Cholesky factor (corresponding to constructor call)
	sparse_d_matrix L_;

	// Number of non-zeros and sparsity pattern for L (set by constructor)
	s_vector                      L_nnz_;
	CppAD::mixed::sparse_mat_info L_pattern_;

	// Forward and reverse values for Alow and L
	CppAD::vector<sparse_d_matrix> f_Alow_, f_L_;
	// -----------------------------------------------------------------
	// functions
	// -----------------------------------------------------------------
	// get the number of non-zeros corresponding to a matrix
	template <class Eigen_sparse_matrix_type>
	s_vector get_nnz(const Eigen_sparse_matrix_type& mat)
	{	typedef typename Eigen_sparse_matrix_type::InnerIterator iterator;
		size_t nc = size_t( mat.cols() );
		s_vector nnz(nc);
		for(size_t j = 0; j < nc; j++)
		{	nnz[j] = 0;
			for(iterator itr(mat, j); itr; ++itr)
			{	nnz[j]++;
			}
		}
		return nnz;
	}
	// -----------------------------------------------------------------
	// use values of elements in Alow to compute new values for L and P
	// (Note that its is expected that P will not change)
	bool factor(
		const sparse_d_matrix& Alow ,
		sparse_d_matrix&       L    ,
		perm_matrix&           P    )
	{	ldlt_obj_.factorize( Alow );
		if( ldlt_obj_.info() != Eigen::Success )
			return false;
		Eigen::VectorXd D2 = sqrt( ldlt_obj_.vectorD().array() ).matrix();
		L                  =  ldlt_obj_.matrixL() * D2.asDiagonal();
		P                  = ldlt_obj_.permutationP();
		return true;
	}
	// ------------------------------------------------------------------
	// get the sparsity pattern corresponding to a matrix
	CppAD::mixed::sparse_mat_info get_pattern(const sparse_d_matrix& mat)
	{	size_t nc = mat.cols();
		CppAD::mixed::sparse_mat_info pattern;
		assert( pattern.row.size() == 0 );
		assert( pattern.col.size() == 0 );
		for(size_t j = 0; j < nc; j++)
		{	int previous = mat.rows();
			for(sparse_d_matrix::InnerIterator itr(mat, j); itr; ++itr)
			{	assert( itr.row() >= 0  && itr.row() < mat.rows());
				assert( size_t(itr.col()) == j );
				pattern.row.push_back( size_t( itr.row() ) );
				pattern.col.push_back(j);
				assert( previous == mat.rows() || previous < itr.row() );
				previous = itr.row();
			}
		}
		return pattern;
	}
	// ------------------------------------------------------------------
	// CppAD forward mode for this operation
	// (Only order zero implemented so far)
	virtual bool forward(
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
	{	typedef typename sparse_d_matrix::InnerIterator iterator;
		assert( p <= q );
		assert( q == 0 );
		//
		size_t n_order = q + 1;
		size_t nx      = Alow_pattern_.row.size();
		size_t ny      = L_pattern_.row.size();
		size_t nc      = Alow_nnz_.size();
		assert( vx.size() == 0 || nx == vx.size() );
		assert( vx.size() == 0 || ny == vy.size() );
		assert( nx * n_order == tx.size() );
		assert( ny * n_order == ty.size() );
		// -------------------------------------------------------------------
		// make sure f_Alow_ and f_L_ are large enough
		assert( f_Alow_.size() == f_L_.size() );
		if( f_Alow_.size() < n_order )
		{	f_Alow_.resize(n_order);
			f_L_.resize(n_order);
			//
			for(size_t k = 0; k < n_order; k++)
			{	f_Alow_[k].resize(nc, nc);
				f_Alow_[k].reserve(Alow_nnz_);
				// indices for possibly non-zero elements in Alow
				for(size_t ell = 0; ell < Alow_pattern_.row.size(); ell++)
				{	size_t i = Alow_pattern_.row[ell];
					size_t j = Alow_pattern_.col[ell];
					f_Alow_[k].insert(i, j) = 1.;
				}
			}
		}
		// -------------------------------------------------------------------
		// unpack tx into f_Alow_
		for(size_t k = 0; k < n_order; k++)
		{	size_t index = 0;
			// unpack Alow values for this order
			for(size_t j = 0; j < nc; j++)
			{	for(iterator itr(f_Alow_[k], j); itr; ++itr)
				{	itr.valueRef() = tx[ index * n_order + k ];
					index++;
				}
			}
		}
		// -------------------------------------------------------------------
		// result value for order k
		ldlt_obj_.factorize( f_Alow_[0] );
		if( ldlt_obj_.info() != Eigen::Success )
		{	ok_ = false;
			return false;
		}
		Eigen::VectorXd D2 = sqrt( ldlt_obj_.vectorD().array() ).matrix();
		f_L_[0]            =  ldlt_obj_.matrixL() * D2.asDiagonal();
		assert( P_.indices() == ldlt_obj_.permutationP().indices() );
		// -------------------------------------------------------------------
		// pack f_L_ into ty
		// -------------------------------------------------------------------
		for(size_t k = 0; k < n_order; k++)
		{	size_t index = 0;
			for(size_t j = 0; j < nc; j++)
			{	for(iterator itr(f_L_[k], j); itr; ++itr)
				{	ty[ index * n_order + k ] = itr.value();
					++index;
				}
			}
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
	// -----------------------------------------------------------------
	// constructor
	atomic_ad_cholesky(const sparse_d_matrix& Alow )
	: CppAD::atomic_base<double> (
		"atomic_ad_cholesky",
		CppAD::atomic_base<double>::set_sparsity_enum
	)
	{	ok_ = true;
		//
		// Step 1: Set Alow_nnz_ and Alow_pattern_
		Alow_nnz_     = get_nnz(Alow);
		Alow_pattern_ = get_pattern(Alow);
		//
		// Step 2: analyze the sparsity pattern
		ldlt_obj_.analyzePattern( Alow );
		//
		// Step 3: Compute the factor L and permutation P for this Alow
		ok_ &= factor(Alow, L_, P_);
		//
		// Step 4: Set L_nnz_ and L_pattern_
		L_nnz_      = get_nnz(L_);
		L_pattern_  = get_pattern(L_);
	}
	// -----------------------------------------------------------------
	// user AD version of atomic Cholesky factorization
	void ad(
		const sparse_ad_matrix& aAlow  ,
		sparse_ad_matrix&       aL     )
	{	// -----------------------------------------------------------
		// packed version of Alow
		size_t nc = Alow_nnz_.size();
		size_t nx = Alow_pattern_.row.size();
		assert( nc == size_t( aAlow.rows() ) );
		assert( nc == size_t( aAlow.cols() ) );
		ad_vector ax( nx );
		size_t index = 0;
		for(size_t j = 0; j < nc; j++)
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
		ad_vector ay( ny );
		(*this)(ax, ay);
		// -------------------------------------------------------------------
		// unpack ay into aL
		aL.resize(nc, nc);
		index = 0;
		for(size_t ell = 0; ell < ny; ell++)
		{	size_t i = L_pattern_.row[ell];
			size_t j = L_pattern_.col[ell];
			aL.insert(i, j) = ay[ index++ ];
		}
		assert( index == ny );
		return;
	}
}; // END_ATOMIC_AD_CHOLESKY

} // END_EMPTY_NAMESPACE

/*
$begin ad_cholesky.cpp$$

$section Sparse Atomic AD Cholesky Factorization: Example and Test$$

$head Problem$$
We are given the function $latex A : \B{R}^3 \rightarrow \B{R}^{3 \times 3}$$
defined by
$latex \[
	A(x) = \left( \begin{array}{ccc}
		x_0 & 0    & x_1 \\
		0   & x_1  & 0   \\
		x_1 & 0    & x_2
	\end{array} \right)
\] $$
The leading princial minors of this matrix are
$latex x_0$$,
$latex x_0 x_1$$,
$latex x_0 x_1 x_2 - x_1 x_1 x_1$$,
If all these minors are positive, the matrix $latex A(x)$$ is
positive definite.

$end
*/

bool ad_cholesky(void)
{	bool ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();
	//
	// create ad_cholesky object
	size_t nc = 3;
	sparse_d_matrix Blow(nc, nc);
	Blow.insert(0,0) = 1.0;
	Blow.insert(2,0) = 1.0;
	Blow.insert(1,1) = 1.0;
	Blow.insert(2,2) = 2.0;
	atomic_ad_cholesky cholesky( Blow );
	//
	// Independent variables
	size_t nx = 3;
	ad_vector ax(nx);
	ax[0] = 2.0;
	ax[1] = 1.0;
	ax[2] = 3.0;
	CppAD::Independent( ax );
	//
	// Lower triangle of symmetric matrix with same sparsity pattern as B
	sparse_ad_matrix aAlow(nc, nc);
	aAlow.insert(0,0) = ax[0];
	aAlow.insert(2,0) = ax[1];
	aAlow.insert(1,1) = ax[1];
	aAlow.insert(2,2) = ax[2];
	//
	// compute the Choleksy factorization of A
	sparse_ad_matrix aL;
	cholesky.ad(aAlow, aL);
	//
	// diagonal of L
	Eigen::Matrix< AD<double> , Dynamic , 1 > D = aL.diagonal();
	//
	// product of diagonal elements of L
	AD<double> p = 1.0;
	for(size_t j = 0; j < nc; j++)
		p *= D[j];
	//
	// calculate the determinant of A using the Cholesky factorization
	ad_vector ay(1);
	ay[0] = p * p;
	//
	// f(x) = det[ A(x) ]
	CppAD::ADFun<double> f(ax, ay);
	//
	// Check the determinant
	AD<double> acheck = ax[0] * ax[1] * ax[2] - ax[1] * ax[1] * ax[1];
	ok &= CppAD::abs( acheck - ay[0] ) < eps;
	//
	// Check that A is equal to  P' * L * L' * P
	sparse_d_matrix L     = admat2mat( aL , cholesky.L_nnz_ );
	sparse_d_matrix LLT   = L * L.transpose();
	sparse_d_matrix PTLLT = cholesky.P_.transpose() * LLT;
	sparse_d_matrix prod  = PTLLT * cholesky.P_.transpose();
	//
	dense_d_matrix temp = prod;
	dense_d_matrix A(3,3);
	A << 2, 0, 1, 0, 1, 0, 1, 0, 3;
	for(size_t i = 0; i < nc; i++)
	{	for(size_t j = 0; j < nc; j++)
		{	ok &= std::fabs( A(i, j) - temp(i, j) ) < eps;
		}
	}
	// -----------------------------------------------------------------------
	ok &= cholesky.ok_;
	return ok;
}
