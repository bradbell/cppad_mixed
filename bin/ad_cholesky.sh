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
# include <iostream>


int main()
{	typedef Eigen::Matrix<
		double, Eigen::Dynamic, Eigen::Dynamic>        dense_matrix;
	typedef Eigen::PermutationMatrix<
		Eigen::Dynamic, Eigen::Dynamic>                perm_matrix;
	typedef Eigen::SparseMatrix<double>                sparse_matrix;
	typedef Eigen::SparseMatrix< CppAD::AD<double> >   ad_sparse_matrix;
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
	// Step 1: determine number of non-zeros in the columns of aAlow
	Eigen::VectorXi nnz(nc);
	for(size_t j = 0; j < nc; j++)
	{	nnz[j] = 0;
		for(ad_sparse_matrix::InnerIterator itr(aAlow, j); itr; ++itr)
		{	nnz[j]++;
			assert( itr.row() >= itr.col() );
		}
	}
	//
	// Step 2: set Alow to a double version of aAlow 
	sparse_matrix Alow(nc, nc);
	Alow.reserve(nnz);
	for(size_t j = 0; j < nc; j++)
	{	nnz[j] = 0;
		for(ad_sparse_matrix::InnerIterator itr(aAlow, j); itr; ++itr)
		{	CppAD::AD<double> av = itr.value();
			av                   = Var2Par( av );
			double            v  = Value( av );
			Alow.insert( itr.row(), itr.col() ) = v;
		}
	}
	//
	// Step 3: compute factorization
	Eigen::SimplicialLDLT<sparse_matrix> ldlt( Alow );
	perm_matrix P     = ldlt.permutationP();
	sparse_matrix L   = ldlt.matrixL();
	Eigen::VectorXd D = ldlt.vectorD();
	// -----------------------------------------------------------------------
	// Test: check that A is equal to  P' * L * D * L' * P 
	sparse_matrix prod = L * D.asDiagonal() * L.transpose();
	//
	// Cannot seem to get this operation to work with sparse matrices
	// prod  = P.transpose() * prod * P;
	// For now, convert to dense matrices where it does work
	dense_matrix temp = prod;
	temp = P.transpose() * temp * P;
	std::cout << "P'*L*D*L'*P =\n" << temp << std::endl;
	//
	// Test for A
	dense_matrix B = Alow;
	dense_matrix A = B + B.transpose();
	for(size_t j = 0; j < nc; j++)
		A(j, j) = A(j, j) / 2.0;
	dense_matrix check(3,3);
	check << 2, 0, 1, 0, 1, 0, 1, 0, 2;
	ok = check == A;
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
	-I$HOME/prefix/eigen/include \
	-I$HOME/prefix/cppad/include \
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
