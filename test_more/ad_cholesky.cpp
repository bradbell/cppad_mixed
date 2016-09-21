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
// define some types of matrices
// -----------------------------------------------------------------
typedef Eigen::Matrix<
	double, Eigen::Dynamic, Eigen::Dynamic>        dense_matrix;
typedef Eigen::PermutationMatrix<
	Eigen::Dynamic, Eigen::Dynamic>                perm_matrix;
typedef Eigen::SparseMatrix<double>                sparse_matrix;
typedef Eigen::SparseMatrix< CppAD::AD<double> >   ad_sparse_matrix;

class atomic_ad_cholesky : public CppAD::atomic_base<double> {
public:
	// -----------------------------------------------------------------
	// data (in order that it is initialized)
	// -----------------------------------------------------------------
	// OK flag
	bool ok_;
	//
	// Number of non-zeros and sparsity pattern for Alow (set by constructor)
	Eigen::VectorXi               Alow_nnz_;
	CppAD::mixed::sparse_mat_info Alow_pattern_;

	// Object used for Cholesky factorization
	// (analyzePattern by constructor only)
	Eigen::SimplicialLDLT<sparse_matrix> ldlt_obj_;

	// Value of the permutation matrix (set by constructor)
	perm_matrix P_;

	// Value of the Cholesky factor (corresponding to constructor call)
	sparse_matrix L_;

	// Number of non-zeros and sparsity pattern for L (set by constructor)
	Eigen::VectorXi               L_nnz_;
	CppAD::mixed::sparse_mat_info L_pattern_;

	// Forward and reverse values for Alow and L
	CppAD::vector<sparse_matrix> f_Alow_, f_L_;
	// -----------------------------------------------------------------
	// functions
	// -----------------------------------------------------------------
	// get the number of non-zeros corresponding to a matrix
	template <class Eigen_sparse_matrix_type>
	Eigen::VectorXi get_nnz(const Eigen_sparse_matrix_type& mat)
	{	typedef typename Eigen_sparse_matrix_type::InnerIterator iterator;
		size_t nc = size_t( mat.cols() );
		Eigen::VectorXi nnz(nc);
		for(size_t j = 0; j < nc; j++)
		{	nnz[j] = 0;
			for(iterator itr(mat, j); itr; ++itr)
			{	nnz[j]++;
			}
		}
		return nnz;
	}
	// -----------------------------------------------------------------
	// convert ad_sparse_matrix -> sparse_matrix
	sparse_matrix amat2dmat(
		const ad_sparse_matrix& amat  ,
		const Eigen::VectorXi&  nnz   )
	{	size_t nr = size_t( amat.rows() );
		size_t nc = size_t( amat.cols() );
		assert( size_t( nnz.size() ) == nc );
		//
		sparse_matrix dmat(nr, nc);
		dmat.reserve(nnz);
		for(size_t j = 0; j < nc; j++)
		{	for(ad_sparse_matrix::InnerIterator itr(amat, j); itr; ++itr)
			{	CppAD::AD<double> av = itr.value();
				av                   = Var2Par( av );
				double            v  = Value( av );
				dmat.insert( itr.row(), itr.col() ) = v;
			}
		}
		return dmat;
	}
	// -----------------------------------------------------------------
	// use values of elements in Alow to compute new values for L and P
	// (Note that its is expected that P will not change)
	bool factor(const sparse_matrix& Alow, sparse_matrix& L, perm_matrix& P)
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
	CppAD::mixed::sparse_mat_info get_pattern(const sparse_matrix& mat)
	{	size_t nc = mat.cols();
		CppAD::mixed::sparse_mat_info pattern;
		assert( pattern.row.size() == 0 );
		assert( pattern.col.size() == 0 );
		for(size_t j = 0; j < nc; j++)
		{	int previous = mat.rows();
			for(sparse_matrix::InnerIterator itr(mat, j); itr; ++itr)
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
	{	typedef typename sparse_matrix::InnerIterator iterator;
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
			{	f_Alow_[k].setZero();
				f_L_[k].setZero();
				//
				f_Alow_[k].resize(nc, nc);
				f_L_[k].resize(nc, nc);
				//
				f_Alow_[k].reserve(Alow_nnz_);
				f_L_[k].reserve(L_nnz_);
				//
				// indices for possibly non-zero elements in Alow
				for(size_t ell = 0; ell < Alow_pattern_.row.size(); ell++)
				{	size_t i = Alow_pattern_.row[ell];
					size_t j = Alow_pattern_.col[ell];
					f_Alow_[k].insert(i, j) = 0.;
				}
				//
				// indices for possibly non-zero elements in L
				for(size_t ell = 0; ell < L_pattern_.row.size(); ell++)
				{	size_t i = L_pattern_.row[ell];
					size_t j = L_pattern_.col[ell];
					f_L_[k].insert(i, j) = 0.;
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
			return false;
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
					ty[ index * n_order + k ] = itr.value();
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
	atomic_ad_cholesky(const sparse_matrix& Alow )
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
}; // END_ATOMIC_AD_CHOLESKY

} // END_EMPTY_NAMESPACE

bool ad_cholesky(void)
{	bool ok = true;
	//
	double eps = 100. * std::numeric_limits<double>::epsilon();
	//
	// AD version of lower triangle of the symmetric positive difinite matrix
	// A = [ 2, 0, 1; 0, 1, 0; 1. 0, 2 ];
	size_t nc = 3;
	sparse_matrix Alow(nc, nc);
	Alow.insert(0,0) = 2.0;
	Alow.insert(1,1) = 1.0;
	Alow.insert(2,0) = 1.0;
	Alow.insert(2,2) = 2.0;
	//
	// object that computes the choleksy factor of a matrix
	atomic_ad_cholesky cholesky( Alow );
	ok &= cholesky.ok_;
	// -----------------------------------------------------------------------
	// Test sparsity pattern for Alow is in column major order
	ok &= cholesky.Alow_pattern_.row.size() == 4;
	ok &= cholesky.Alow_pattern_.col.size() == 4;
	//
	ok &= cholesky.Alow_pattern_.row[0]     == 0;
	ok &= cholesky.Alow_pattern_.col[0]     == 0;
	//
	ok &= cholesky.Alow_pattern_.row[1]     == 2;
	ok &= cholesky.Alow_pattern_.col[1]     == 0;
	//
	ok &= cholesky.Alow_pattern_.row[2]     == 1;
	ok &= cholesky.Alow_pattern_.col[2]     == 1;
	//
	ok &= cholesky.Alow_pattern_.row[3]     == 2;
	ok &= cholesky.Alow_pattern_.col[3]     == 2;
	// -----------------------------------------------------------------------
	// Test: check that A is equal to  P' * L * L' * P
	sparse_matrix LLT   = cholesky.L_ * cholesky.L_.transpose();
	sparse_matrix PTLLT = cholesky.P_.transpose() * LLT;
	sparse_matrix prod   = PTLLT * cholesky.P_.transpose();
	//
	dense_matrix temp = prod;
	dense_matrix A(3,3);
	A << 2, 0, 1, 0, 1, 0, 1, 0, 2;
	for(size_t i = 0; i < nc; i++)
	{	for(size_t j = 0; j < nc; j++)
			ok &= std::fabs( A(i, j) - temp(i, j) ) < eps;
	}
	// -----------------------------------------------------------------------
	return ok;
}
