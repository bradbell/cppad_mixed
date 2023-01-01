/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin ldlt_eigen.cpp$$
$spell
	rcv
	rc
	nrow
	init
	CppAD
	eigen ldlt_obj
	logdet
	sim_cov
$$

$nospell
$bold This is cppad_mixed--20220519 documentation:$$ Here is a link to its
$href%https://cppad-mixed.readthedocs.io%current documentation%$$.
$$
$section Example Using Eigen LDLT Factorization$$

$head Problem Description$$
We define the lower triangular matrix
$latex \[
	L =
	\left( \begin{array}{ccc}
		1 & 0 & 0 \\
		2 & 1 & 0 \\
		3 & 2 & 1
	\end{array} \right)
	\W{,}
	D =
	\left( \begin{array}{ccc}
		3 & 0 & 0 \\
		0 & 2 & 0 \\
		0 & 0 & 1
	\end{array} \right)
	\W{,}
	L^\R{T} =
	\left( \begin{array}{ccc}
		1 & 2 & 3 \\
		0 & 1 & 2 \\
		0 & 0 & 1
	\end{array} \right)
\] $$
and the positive definite matrix
$latex \[
	H = L D L^\R{T} =
	\left( \begin{array}{ccc}
		3 & 6  & 9 \\
		6 & 14 & 22 \\
		9 & 22 & 36
	\end{array} \right)
\] $$
The inverse of $latex H$$ is given by
$latex \[
	H^{-1} =
	\frac{1}{6}
	\left( \begin{array}{ccc}
		20  & -18  & 6   \\
		-18 & 27   & -12 \\
		6   & -12  & 6
	\end{array} \right)
\] $$
which can be checked by multiplying by $latex H H^{-1}$$.

$head constructor$$
See the following code below:
$codep
	CppAD::mixed::ldlt_eigen ldlt_obj(nrow);
$$

$head init$$
See the following under
$cref/Source Code/ldlt_eigen.cpp/Source Code/$$ below:
$codep
	ldlt_obj.init( H_rcv.pat() );
$$

$head pattern$$
See the following under
$cref/Source Code/ldlt_eigen.cpp/Source Code/$$ below:
$codep
	H_rc = ldlt_obj.pattern();
$$

$head update$$
See the following under Source Code below:
$codep
	ldlt_obj.update( H_rcv );
$$

$head logdet$$
See the following under Source Code below:
$codep
	logdet_H = ldlt_obj.logdet(negative);
$$

$head solve_H$$
See the following under Source Code below:
$codep
	ldlt_obj.solve_H(row, val_in, val_out);
$$

$head split$$
See the following under Source Code below:
$codep
	ok &= ldlt_obj.split(L, D, P)
$$

$head solve_LDLT$$
See the following under Source Code below:
$codep
	x = ldlt_eigen<double>::solve_LDLT(L, D, P, b);
$$


$head sim_cov$$
See the following under Source Code below:
$codep
	ok &= ldlt_obj.sim_cov(w, v)
$$

$head Source Code$$
$code
$srcthisfile%5%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/mixed/ldlt_eigen.hpp>
# include <limits>
# include <cmath>
# include <cassert>

bool ldlt_eigen_xam(void)
{	bool ok    = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();

	double H_inv[] = {
		20.0,  -18.0,  6.0,
		-18.0,  27.0, -12.0,
		6.0,   -12.0,  6.0
	};
	for(size_t i = 0; i < sizeof(H_inv)/sizeof(H_inv[0]); i++)
		H_inv[i] /= 6.0;

	// create eigen object
	size_t nrow = 3;    // number of rows in H
	size_t ncol = nrow; // number of columns in H
	CppAD::mixed::ldlt_eigen<double> ldlt_obj(nrow);
	assert( nrow * ncol == sizeof(H_inv) / sizeof(H_inv[0]) );

	// sparsity pattern for the lower triangle of H (dense)
	size_t nnz = 6;
	CppAD::mixed::sparse_rc H_rc(nrow, ncol, nnz);
	{	size_t k = 0;
		H_rc.set(k++, 0, 0);
		H_rc.set(k++, 1, 0);
		H_rc.set(k++, 2, 0);
		H_rc.set(k++, 1, 1);
		H_rc.set(k++, 2, 1);
		H_rc.set(k++, 2, 2);
	}

	// values in lower triangle of H
	CppAD::mixed::d_sparse_rcv H_rcv( H_rc );
	{	size_t k = 0;
		H_rcv.set(k++,  3.0);  // H_0,0 =  3.0
		H_rcv.set(k++,  6.0);  // H_1,0 =  6.0
		H_rcv.set(k++,  9.0);  // H_2,0 =  9.0
		H_rcv.set(k++, 14.0);  // H_1,1 = 14.0
		H_rcv.set(k++, 22.0);  // H_2,1 = 22.0
		H_rcv.set(k++, 36.0);  // H_2,2 = 36.0
	}
	//
	// initialize the matrix using only the sparsity pattern
	ldlt_obj.init( H_rcv.pat() );

	// check the pattern function
	const CppAD::mixed::sparse_rc& pattern( ldlt_obj.pattern() );
	ok &= pattern.nnz() == nnz;
	for(size_t k = 0; k < nnz; ++k)
	{	ok &= pattern.row()[k] == H_rc.row()[k];
		ok &= pattern.col()[k] == H_rc.col()[k];
	}

	// factor the matrix using the values
	ldlt_obj.update( H_rcv );

	// compute log of determinant of H
	size_t negative;
	double logdet_H = ldlt_obj.logdet(negative);
	ok &= negative == 0;

	// check its value
	ok &= std::fabs( logdet_H / std::log(6.0) - 1.0 ) <= eps;

	// test solve_H
	CppAD::vector<size_t> row(3);
	CppAD::vector<double> val_in(3), val_out(3);
	for(size_t j = 0; j < ncol; j++)
	{	// solve for the j-th column of the inverse matrix
		for(size_t k = 0; k < 3; k++)
		{	row[k] = k;
			if( row[k] == j )
				val_in[k] = 1.0;
			else
				val_in[k] = 0.0;
		}
		ldlt_obj.solve_H(row, val_in, val_out);
		//
		for(size_t k = 0; k < row.size(); k++)
		{	size_t i       = row[k];
			double check_i = H_inv[ i * nrow + j ];
			ok &= std::fabs( val_out[k] - check_i ) <= eps;
		}
	}

	// test split
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigen_dense;
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor> eigen_sparse;
	typedef Eigen::Matrix<double, Eigen::Dynamic, 1>     eigen_vector;
	typedef Eigen::PermutationMatrix<Eigen::Dynamic>     eigen_perm;
	eigen_sparse L;
	eigen_vector D;
	eigen_perm   P;
	ldlt_obj.split(L, D, P);
	eigen_dense denL = L;
	eigen_dense denP = P;
	for(size_t i = 0; i < nrow; ++i)
	{	// D is as described above
		ok &= std::fabs( D(i) - double(3 - i) ) <= eps;
		for(size_t j = 0; j < ncol; ++j)
		{	// L is as described above
			double	check = std::max(0.0, double(i) + 1.0 - double(j));
			ok &= std::fabs( denL(i, j) - check ) <= eps;
			// P is identity
			check = double( i == j );
			ok &= std::fabs( denP(i, j) - check ) <= eps;
		}
	}

	// test solve_LDLT
	eigen_vector x(nrow), b(nrow);
	for(size_t j = 0; j < ncol; j++)
	{	// solve for the j-th column of the inverse matrix
		for(size_t i = 0; i < nrow; i++)
		{	if( i == j )
				b[i] = 1.0;
			else
				b[i] = 0.0;
		}
		x = CppAD::mixed::ldlt_eigen<double>::solve_LDLT(L, D, P, b);
		//
		for(size_t i = 0; i < nrow; i++)
		{	double check_i = H_inv[ i * nrow + j ];
			ok &= std::fabs( x[i] - check_i ) <= eps;
		}
	}

	// test sim_cov
	CppAD::vector<double> w(3), v(3), c(3);
	for(size_t i = 0; i < 3; i++)
		w[i] = double(2 * i + 1);
	ok &= ldlt_obj.sim_cov(w, v);
	// check that w =  sqrt(D) * L^T P * v = sqrt(D) * L^T * v
	for(size_t i = 0; i < 3; i++)
	{	double check = 0.0;
		for(size_t j = i; j < 3; j++)
			check += denL(j, i) * v[j];
		check *= std::sqrt( D[i] );
		ok &= std::fabs( w[i] - check ) <= eps;
	}

	return ok;
}
// END C++
