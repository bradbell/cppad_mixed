/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/ldlt_cholmod.hpp>
# include <limits>
# include <cmath>
# include <cassert>

namespace { // BEGIN_EMPTY_NAMESPACE
/*
begin ldlt_cholmod_xam.cpp
$spell
	nrow
	init
	Cholesky
	CppAD
	cholmod ldlt_obj
	logdet
$$

$nospell
$bold This is cppad_mixed--20220519 documentation:$$ Here is a link to its
$href%https://cppad-mixed.readthedocs.io%current documentation%$$.
$$
$section Example Using Cholmod Cholesky Factorization$$

$head Problem Description$$
Solve for $latex x$$ in the equation $latex H x = b$$ where $latex H$$
is defined below and $latex b$$ is a column of the identity matrix.
Hence the solution $latex x$$ is the corresponding column of
$latex H^{-1}$$.
$latex \[
	G = \left( \begin{array}{ccc}
		5 & 4 & 2 \\
		4 & 5 & 1 \\
		2 & 1 & 5
	\end{array} \right)
	\W{,}
	H = \left( \begin{array}{cc}
		G & 0 \\
		0 & G
	\end{array} \right)
\] $$
We use $latex G^k$$ to denote the upper-left $latex k \times k$$
principal minor. The determinant of its principal minors are:
$latex \[
	\det \left( G^1 \right) = 5  \W{,}
	\det \left( G^2 \right) = 9  \W{,}
	\det \left( G^3 \right) = 36
\] $$
Hence, $latex G$$ and $latex H$$ are positive definite.
In addition
$latex \[
	G^{-1} = \frac{1}{36}
	\left( \begin{array}{ccc}
		24  & -18 & -6 \\
		-18 & 21  & 3 \\
		-6  & 3   & 9
	\end{array} \right)
	\W{,}
	H^{-1} = \left( \begin{array}{cc}
		G^{-1} & 0 \\
		0 & G^{-1}
	\end{array} \right)
\] $$
which can be checked by multiplying by $latex G G^{-1}$$.

$head constructor$$
See the following code below:
$codep
	CppAD::mixed::ldlt_cholmod ldlt_obj(nrow);
$$

$head init$$
See the following under
$cref/Source Code/ldlt_cholmod.cpp/Source Code/$$ below:
$codep
	ldlt_obj.init( H_rcv.pat() );
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

$head Source Code$$
$code
$srcthisfile%5%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++

bool ldlt_cholmod_1(void)
{	bool ok    = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();

	double H_inv[] = {
		 24.0, -18.0, -6.0,   0.0,   0.0,  0.0,
		-18.0,  21.0,  3.0,   0.0,   0.0,  0.0,
		 -6.0,   3.0,  9.0,   0.0,   0.0,  0.0,
		  0.0,   0.0,  0.0,  24.0, -18.0, -6.0,
		  0.0,   0.0,  0.0, -18.0,  21.0,  3.0,
		  0.0,   0.0,  0.0,  -6.0,   3.0,  9.0
	};
	for(size_t i = 0; i < sizeof(H_inv)/sizeof(H_inv[0]); i++)
		H_inv[i] /= 36.;

	// create cholmod object
	size_t nrow = 6;    // number of rows in H
	size_t ncol = nrow; // number of columns in H
	CppAD::mixed::ldlt_cholmod ldlt_obj(nrow);
	assert( nrow * ncol == sizeof(H_inv) / sizeof(H_inv[0]) );

	// sparsity pattern corresponding to lower triangle of H
	// in column major order
	size_t nnz = 12;
	CppAD::mixed::sparse_rc H_rc(nrow, ncol, nnz);
	{	size_t k = 0;
		// upper left block
		H_rc.set(k++, 0, 0);
		H_rc.set(k++, 1, 0);
		H_rc.set(k++, 2, 0);
		H_rc.set(k++, 1, 1);
		H_rc.set(k++, 2, 1);
		H_rc.set(k++, 2, 2);
		// lower right block
		H_rc.set(k++, 3, 3);
		H_rc.set(k++, 4, 3);
		H_rc.set(k++, 5, 3);
		H_rc.set(k++, 4, 4);
		H_rc.set(k++, 5, 4);
		H_rc.set(k++, 5, 5);
	}
	// values in lower triangle of H
	CppAD::mixed::d_sparse_rcv H_rcv( H_rc );
	{	size_t k = 0;
		H_rcv.set(k++, 5.0); // H_0,0  = 5.0
		H_rcv.set(k++, 4.0); // H_1,0  = 4.0
		H_rcv.set(k++, 2.0); // H_2,0  = 2.0
		H_rcv.set(k++, 5.0); // H_1,1  = 5.0
		H_rcv.set(k++, 1.0); // H_2,1  = 1.0
		H_rcv.set(k++, 5.0); // H_2,2  = 5.0
		H_rcv.set(k++, 5.0); // H_3,3  = 5.0
		H_rcv.set(k++, 4.0); // H_4,3  = 4.0
		H_rcv.set(k++, 2.0); // H_5,3  = 2.0
		H_rcv.set(k++, 5.0); // H_4,4  = 5.0
		H_rcv.set(k++, 1.0); // H_5,4  = 1.0
		H_rcv.set(k++, 5.0); // H_5,5  = 5.0
	}

	// initialize the matrix using only the sparsity pattern
	ldlt_obj.init( H_rcv.pat() );

	// factor the matrix using the values
	ldlt_obj.update( H_rcv );

	// compute log of determinant of H
	size_t negative;
	double logdet_H = ldlt_obj.logdet(negative);
	ok &= negative == 0;

	// check its value
	ok &= std::fabs( logdet_H / (2.0 * std::log(36.0)) - 1.0 ) <= eps;

	// test solve
	{
		CppAD::vector<size_t> row(3);
		CppAD::vector<double> val_in(3), val_out(3);
		for(size_t j = 0; j < ncol; j++)
		{	// solve for the j-th column of the inverse matrix
			for(size_t k = 0; k < 3; k++)
			{	if( j < 3 )
					row[k] = k;
				else
					row[k] = k + 3;
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
				ok &= std::fabs( val_out[k] / check_i - 1.0 ) <= eps;
			}
		}
	}
	// test inv on subblock
	{
		size_t K = 9;
		CppAD::vector<size_t> row_in(K), col_in(K);
		CppAD::vector<double> val_out(K);
		{	size_t k = 0;
			// use row major ordering to test sorting
			for(size_t i = 0; i < 3; i++)
			{	for(size_t j = 0; j < 3; j++)
				{	row_in[k] = 3 + i;
					col_in[k] = 3 + j;
					k++;
				}
			}
		}
		ldlt_obj.inv(row_in, col_in, val_out);
		for(size_t k = 0; k < K; k++)
		{	double check = H_inv[ row_in[k] * nrow + col_in[k] ];
			ok          &= std::fabs( val_out[k] - check ) <= eps;
		}
	}
	return ok;
}
// ----------------------------------------------------------------------------
bool ldlt_cholmod_2(void)
{	bool ok    = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();
	//
	// determinant of H is 10
# ifndef NDEBUG
	double H[] = {
		1.0, 0.0, 0.0,  1.0,
		0.0, 1.0, 0.0,  0.0,
		0.0, 0.0, 1.0,  0.0,
		1.0, 0.0, 0.0, 11.0
	};
# endif
	// The inverse of H is
	double H_inv[] = {
		1.1, 0.0, 0.0, -0.1,
		0.0, 1.0, 0.0,  0.0,
		0.0, 0.0, 1.0,  0.0,
		-0.1, 0.0, 0.0, 0.1
	};
	//
	// create cholmod object
	size_t nrow = 4;    // number of rows in H
	size_t ncol = nrow; // number of columns in H
	CppAD::mixed::ldlt_cholmod ldlt_obj(nrow);
	assert( nrow * ncol == sizeof(H) / sizeof(H[0]) );
	//
	// sparsity pattern corresponding to lower triangle of H
	// in column major order
	size_t nnz = 5;
	CppAD::mixed::sparse_rc H_rc(nrow, ncol, nnz);
	{	size_t k = 0;
		// lower triangle
		H_rc.set(k++, 0, 0);
		H_rc.set(k++, 3, 0);
		H_rc.set(k++, 1, 1);
		H_rc.set(k++, 2, 2);
		H_rc.set(k++, 3, 3);
	}
	// values in lower triangle of H
	CppAD::mixed::d_sparse_rcv H_rcv( H_rc );
	{	size_t k = 0;
		H_rcv.set(k++,  1.0); // H_0,0  =  1.0
		H_rcv.set(k++,  1.0); // H_3,0  =  1.0
		H_rcv.set(k++,  1.0); // H_1,1  =  1.0
		H_rcv.set(k++,  1.0); // H_2,2  =  1.0
		H_rcv.set(k++, 11.0); // H_2,2  = 11.0
	}
	//
	// initialize the matrix using only the sparsity pattern
	ldlt_obj.init( H_rcv.pat() );
	//
	// factor the matrix using the values
	ldlt_obj.update( H_rcv );
	//
	// compute log of determinant of H
	size_t negative;
	double logdet_H = ldlt_obj.logdet(negative);
	ok &= negative == 0;
	//
	// check its value
	ok &= std::fabs( logdet_H / std::log(10.0) - 1.0 ) <= eps;
	//
	// use solve_H to solve for the entie inverse matrix one column at a time
	{	CppAD::vector<size_t> row(4);
		for(size_t k = 0; k < 4; ++k)
			row[k] = k;
		CppAD::vector<double> val_in(4), val_out(4);
		for(size_t j = 0; j < ncol; j++)
		{	// solve for the j-th column of the inverse matrix
			for(size_t k = 0; k < 4; k++)
			{	val_in[k] = 0.0;
				if( k == j )
					val_in[k] = 1.0;
			}
			ldlt_obj.solve_H(row, val_in, val_out);
			//
			for(size_t k = 0; k < row.size(); k++)
			{	size_t i       = row[k];
				double check_i = H_inv[ i * nrow + j ];
				ok &= std::fabs( val_out[k] - check_i ) <= eps;
			}
		}
	}
	{
		CppAD::vector<size_t> row_in(nnz), col_in(nnz);
		CppAD::vector<double> val_out(nnz);
		{	size_t k = 0;
			for(size_t i = 0; i < 4; i++)
			{	row_in[k] = i;
				col_in[k] = i;
				k++;
			}
			row_in[4] = 3;
			col_in[4] = 0;
		}
		ldlt_obj.inv(row_in, col_in, val_out);
		for(size_t k = 0; k < nnz; k++)
		{	double check = H_inv[ row_in[k] * nrow + col_in[k] ];
			ok          &= std::fabs( val_out[k] - check ) <= eps;
		}
	}
	return ok;
}

} // END_EMPTY_NAMESPACE

bool ldlt_cholmod(void)
{	bool ok = true;
	ok     &= ldlt_cholmod_1();
	ok     &= ldlt_cholmod_2();
	return ok;
}
// END C++
