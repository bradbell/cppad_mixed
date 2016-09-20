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
	// define some types of matrices
	typedef Eigen::Matrix<
		double, Eigen::Dynamic, Eigen::Dynamic>        dense_matrix;
	typedef Eigen::PermutationMatrix<
		Eigen::Dynamic, Eigen::Dynamic>                perm_matrix;
	typedef Eigen::SparseMatrix<double>                sparse_matrix;
	typedef Eigen::SparseMatrix< CppAD::AD<double> >   ad_sparse_matrix;

	// Sparsity pattern for Alow and L in column major order
	CppAD::mixed::sparse_mat_info Alow_pattern, L_pattern;

	// Object used for Cholesky factorization:
	static Eigen::SimplicialLDLT<sparse_matrix> ldlt_obj;

	bool factor(const sparse_matrix& Alow, sparse_matrix& L, perm_matrix& P)
	{	ldlt_obj.factorize( Alow );
		if( ldlt_obj.info() != Eigen::Success )
			return false;
		Eigen::VectorXd D2 = sqrt( ldlt_obj.vectorD().array() ).matrix();
		P                  = ldlt_obj.permutationP();
		L                  =  ldlt_obj.matrixL() * D2.asDiagonal();
		return true;
	}

	// convert ad_sparse_matrix -> sparse_matrix
	sparse_matrix amat2dmat(const ad_sparse_matrix& amat)
	{	size_t nr = size_t( amat.rows() );
		size_t nc = size_t( amat.cols() );
		//
		// determine number of non-zeros in the columns of amat
		Eigen::VectorXi nnz(nc);
		for(size_t j = 0; j < nc; j++)
		{	nnz[j] = 0;
			for(ad_sparse_matrix::InnerIterator itr(amat, j); itr; ++itr)
			{	nnz[j]++;
			}
		}
		// set dmat to a double version of amat
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

	// set the sparsity pattern corresponding to a matrix
	void set_pattern(
		const sparse_matrix&            mat     ,
		CppAD::mixed::sparse_mat_info&  pattern )
	{	size_t nc = mat.cols();
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
	ad_sparse_matrix aAlow(nc,nc);
	aAlow.insert(0,0) = 2.0;
	aAlow.insert(1,1) = 1.0;
	aAlow.insert(2,0) = 1.0;
	aAlow.insert(2,2) = 2.0;
	//
	// Step 1: Convert Alow from AD to double
	sparse_matrix Alow = amat2dmat(aAlow);
	//
	// Step 2: analyze the sparsity pattern
	ldlt_obj.analyzePattern( Alow );
	//
	// Step 3: Compute the factor L and permutation P for this Alow
	sparse_matrix L;
	perm_matrix   P;
	ok &= factor(Alow, L, P);
	//
	// Step 4: Set the sparsity pattern corresponding to Alow and L
	set_pattern(Alow, Alow_pattern);
	set_pattern(L,    L_pattern);
	// -----------------------------------------------------------------------
	// Test sparsity pattern for Alow is in column major order
	ok &= Alow_pattern.row.size() == 4;
	ok &= Alow_pattern.col.size() == 4;
	//
	ok &= Alow_pattern.row[0]     == 0;
	ok &= Alow_pattern.col[0]     == 0;
	//
	ok &= Alow_pattern.row[1]     == 2;
	ok &= Alow_pattern.col[1]     == 0;
	//
	ok &= Alow_pattern.row[2]     == 1;
	ok &= Alow_pattern.col[2]     == 1;
	//
	ok &= Alow_pattern.row[3]     == 2;
	ok &= Alow_pattern.col[3]     == 2;
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
