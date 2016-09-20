#! /bin/bash -e
# $Id$
# ----------------------------------------------------------------------------
# Test methods needed for AD Atomic Cholesky Factorzation
# ----------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-16 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$0" != "bin/ad_cholesky.sh" ]
then
	echo "bin/ad_cholesky.sh: must be executed from its parent directory"
	exit 1
fi
# -----------------------------------------------------------------------------
set -v
if [ ! -e build/ad_cholesky ]
then
	mkdir -p build/ad_cholesky
fi
cd build/ad_cholesky
# ---------------------------------------------------------------------------
set +v
echo 'create build/ad_choleksy/ad_cholesky.cpp'
cat << EOF > ad_cholesky.cpp
# include <Eigen/SparseCore>
# include <Eigen/SparseCholesky>
# include <cppad/cppad.hpp>
# include <cppad/mixed/sparse_mat_info.hpp>
# include <iostream>

namespace {
	// -----------------------------------------------------------------
	// define some types of matrices
	// -----------------------------------------------------------------
	typedef Eigen::Matrix<
		double, Eigen::Dynamic, Eigen::Dynamic>        dense_matrix;
	typedef Eigen::PermutationMatrix<
		Eigen::Dynamic, Eigen::Dynamic>                perm_matrix;
	typedef Eigen::SparseMatrix<double>                sparse_matrix;
	typedef Eigen::SparseMatrix< CppAD::AD<double> >   ad_sparse_matrix;

	// -----------------------------------------------------------------
	// data (in order that it is initialized)
	// -----------------------------------------------------------------

	// Number of non-zeros in each column of Alow (set once)
	Eigen::VectorXi Alow_nnz_;

	// Object used for Cholesky factorization (analyzePattern once)
	static Eigen::SimplicialLDLT<sparse_matrix> ldlt_obj_;

	// Number of non-zeros in each column of L (set once)
	Eigen::VectorXi L_nnz_;

	// Sparsity pattern for Alow and L in column major order (set once)
	CppAD::mixed::sparse_mat_info Alow_pattern_, L_pattern_;

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

	bool factor(const sparse_matrix& Alow, sparse_matrix& L, perm_matrix& P)
	{	ldlt_obj_.factorize( Alow );
		if( ldlt_obj_.info() != Eigen::Success )
			return false;
		Eigen::VectorXd D2 = sqrt( ldlt_obj_.vectorD().array() ).matrix();
		P                  = ldlt_obj_.permutationP();
		L                  =  ldlt_obj_.matrixL() * D2.asDiagonal();
		return true;
	}

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
}

int main()
{
	//
	double eps = 100. * std::numeric_limits<double>::epsilon();
	//
	// Initialize erorr flag
	bool ok = true;
	//
	// AD version of lower triangle of the symmetric positive difinite matrix
	// A = [ 2, 0, 1; 0, 1, 0; 1. 0, 2 ];
	size_t nc = 3;
	ad_sparse_matrix aAlow(nc, nc);
	aAlow.insert(0,0) = 2.0;
	aAlow.insert(1,1) = 1.0;
	aAlow.insert(2,0) = 1.0;
	aAlow.insert(2,2) = 2.0;
	//
	// Step 1: Set Alow_nnz_
	Alow_nnz_ = get_nnz(aAlow);
	//
	// Step 2: Convert Alow from AD to double
	sparse_matrix Alow = amat2dmat(aAlow, Alow_nnz_);
	//
	// Step 3: analyze the sparsity pattern
	ldlt_obj_.analyzePattern( Alow );
	//
	// Step 4: Compute the factor L and permutation P for this Alow
	sparse_matrix L;
	perm_matrix   P;
	ok &= factor(Alow, L, P);
	//
	// Step 5: Set L_nnz_
	L_nnz_ = get_nnz(L);
	//
	// Step 6: Set the sparsity pattern corresponding to Alow and L
	Alow_pattern_ = get_pattern(Alow);
	L_pattern_    = get_pattern(L);
	// -----------------------------------------------------------------------
	// Test sparsity pattern for Alow is in column major order
	ok &= Alow_pattern_.row.size() == 4;
	ok &= Alow_pattern_.col.size() == 4;
	//
	ok &= Alow_pattern_.row[0]     == 0;
	ok &= Alow_pattern_.col[0]     == 0;
	//
	ok &= Alow_pattern_.row[1]     == 2;
	ok &= Alow_pattern_.col[1]     == 0;
	//
	ok &= Alow_pattern_.row[2]     == 1;
	ok &= Alow_pattern_.col[2]     == 1;
	//
	ok &= Alow_pattern_.row[3]     == 2;
	ok &= Alow_pattern_.col[3]     == 2;
	// -----------------------------------------------------------------------
	// Test: check that A is equal to  P' * L * L' * P
	sparse_matrix LLT   = L * L.transpose();
	sparse_matrix PTLLT = P.transpose() * LLT;
	sparse_matrix prod   = PTLLT * P.transpose();
	//
	dense_matrix temp = prod;
	dense_matrix A(3,3);
	A << 2, 0, 1, 0, 1, 0, 1, 0, 2;
	for(size_t i = 0; i < nc; i++)
	{	for(size_t j = 0; j < nc; j++)
			ok &= std::fabs( A(i, j) - temp(i, j) ) < eps;
	}
	// -----------------------------------------------------------------------
	//
	if( ok )
		return 0;
	return 1;
}
EOF
set -v
# ----------------------------------------------------------------------------
g++ \
	-I$HOME/prefix/cppad_mixed/include \
	-I$HOME/prefix/cppad_mixed/eigen/include \
	ad_cholesky.cpp -o ad_cholesky
set +v
if ./ad_cholesky
then
	echo 'ad_cholesky: OK'
else
	echo 'ad_choleksy: Error'
fi
# ----------------------------------------------------------------------------
exit 0
