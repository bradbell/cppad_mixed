// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include "../sparse_ad_cholesky.hpp"

/*
$begin sparse_ad_chol_eq.cpp$$
$spell
	Cholesky
$$

$section Using Sparse AD Cholesky To Solve Equations: Example and Test$$

$tabsize 4$$

$head Source$$
$srcthisfile%4%// BEGIN C++%// END C++%1%$$
$end
*/
// BEGIN C++
bool sparse_ad_chol_eq(void)
{	using CppAD::AD;
	using Eigen::ColMajor;
	using Eigen::Lower;
	using Eigen::Upper;
	typedef Eigen::Matrix< AD<double>, Eigen::Dynamic, 1>  dense_ad_vector;
	//
	bool ok        = true;
	AD<double> eps = 100. * std::numeric_limits<double>::epsilon();
	// --------------------------------------------------------------------
	// create sparse_ad_cholesky object
	int nc = 3;
	Eigen::SparseMatrix< AD<double>, Eigen::ColMajor> ad_Blow(nc, nc);
	ad_Blow.insert(0,0) = 1.0; //     [ 1.0   0.0    0.5 ]
	ad_Blow.insert(2,0) = 0.5; // B = [ 0.0   0.5    0.0 ]
	ad_Blow.insert(1,1) = 0.5; //     [ 0.5   0.0    2.0 ]
	ad_Blow.insert(2,2) = 2.0;
	CppAD::mixed::sparse_ad_cholesky cholesky;
	cholesky.initialize( ad_Blow );
	//
	// Permutation matrix corresponding to this cholesky
	const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>& P =
		cholesky.permutation();
	//
	// Lower triangle of symmetric matrix with same sparsity pattern as B
	Eigen::SparseMatrix< AD<double>, Eigen::ColMajor> Alow(nc, nc);
	Alow.insert(0,0) = 2.0; //     [ 2.0   0.0   0.5  ]
	Alow.insert(2,0) = 0.5; // A = [ 0.0   0.5   0.0  ]
	Alow.insert(1,1) = 0.5; //     [ 0.5   0.0   1.0  ]
	Alow.insert(2,2) = 1.0;
	//
	// compute the Choleksy factorization of A
	Eigen::SparseMatrix< AD<double>, Eigen::ColMajor> L;
	cholesky.eval(Alow, L);
	Eigen::SparseMatrix< AD<double>, Eigen::ColMajor> U = L.transpose();
	ok &= L.rows() == nc;
	ok &= L.cols() == nc;
	//
	// right hand side for equation
	dense_ad_vector b(3);
	b[0] = 1.0;
	b[1] = 2.0;
	b[2] = 3.0;
	//
	// solve the equation: A * x  = b
	//        P * A * P^T * P * x = P * b
	//            L * L^T * P * x = P * b
	dense_ad_vector x, tmp;
	tmp = P * b;
	tmp = L.triangularView<Lower>().solve(tmp);
	tmp = U.triangularView<Upper>().solve(tmp);
	x   = P.transpose() * tmp;
	//
	// dense version of matrix  A
	Eigen::Matrix< AD<double>, 3, 3> A = Alow;
	for(size_t i = 0; i < 3; i++)
	{	for(size_t j = i+1; j < 3; j++)
		{	A(i, j) = A(j, i);
		}
	}
	//
	// compter A * x with b
	dense_ad_vector check = A * x;
	for(size_t i = 0; i < 3; i++)
		ok &= CppAD::NearEqual(check[i], b[i], eps, eps);
	// -----------------------------------------------------------------------
	return ok;
}
// END C++
