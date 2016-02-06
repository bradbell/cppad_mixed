// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin cholmod_xam.cpp$$
$spell
	nrow
	init
	Cholesky
	CppAD
	cholmod cholmod_obj
	logdet
$$

$section Example Using Cholmod Cholesky Factorization$$

$head Problem Description$$
Solve for $latex x$$ in the equation $latex A x = b$$ where $latex A$$
is defined below and $latex b$$ is a column of the identity matrix.
Hence the solution $latex x$$ is the corresponding column of
$latex A^{-1}$$.
$latex \[
	G = \left( \begin{array}{ccc}
		5 & 4 & 2 \\
		4 & 5 & 1 \\
		2 & 1 & 5
	\end{array} \right)
	\W{,}
	A = \left( \begin{array}{cc}
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
Hence, $latex G$$ and $latex A$$ are positive definite.
In addition
$latex \[
	G^{-1} = \frac{1}{36}
	\left( \begin{array}{ccc}
		24  & -18 & -6 \\
		-18 & 21  & 3 \\
		-6  & 3   & 9
	\end{array} \right)
	\W{,}
	A^{-1} = \left( \begin{array}{cc}
		G^{-1} & 0 \\
		0 & G^{-1}
	\end{array} \right)
\] $$
which can be checked by multiplying by $latex G G^{-1}$$.

$head constructor$$
See the following code below:
$codep
	CppAD::mixed::cholmod cholmod_obj(nrow);
$$

$head init$$
See the following under
$cref/Source Code/cholmod_xam.cpp/Source Code/$$ below:
$codep
	cholmod_obj.init(A_info);
$$

$head update$$
See the following under Source Code below:
$codep
	cholmod_obj.update(A_info);
$$

$head logdet$$
See the following under Source Code below:
$codep
	cholmod_obj.logdet(A_info);
$$

$head solve$$
See the following under Source Code below:
$codep
	cholmod_obj.solve(row, val_in, val_out);
$$

$head Source Code$$
$code
$srcfile%example/private/cholmod_xam.cpp%5%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/mixed/cholmod.hpp>
# include <limits>
# include <cmath>
# include <cassert>

bool cholmod_xam(void)
{	bool ok    = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();

	double A_inv[] = {
		 24.0, -18.0, -6.0,   0.0,   0.0,  0.0,
		-18.0,  21.0,  3.0,   0.0,   0.0,  0.0,
		 -6.0,   3.0,  9.0,   0.0,   0.0,  0.0,
		  0.0,   0.0,  0.0,  24.0, -18.0, -6.0,
		  0.0,   0.0,  0.0, -18.0,  21.0,  3.0,
		  0.0,   0.0,  0.0,  -6.0,   3.0,  9.0
	};
	for(size_t i = 0; i < sizeof(A_inv)/sizeof(A_inv[0]); i++)
		A_inv[i] /= 36.;

	// create cholmod object
	size_t nrow = 6;    // number of rows in A
	size_t ncol = nrow; // number of columns in A
	CppAD::mixed::cholmod cholmod_obj(nrow);
	assert( nrow * ncol == sizeof(A_inv) / sizeof(A_inv[0]) );

	// create a sparse matrix representation of the lower triangular of A
	CppAD::mixed::sparse_mat_info A_info;
	A_info.resize(12);
	// A_0,0  = 5.0
	A_info.row[0]  = 0; A_info.col[0]  = 0; A_info.val[0]  = 5.0;
	// A_1,0  = 4.0
	A_info.row[1]  = 1; A_info.col[1]  = 0; A_info.val[1]  = 4.0;
	// A_2,0  = 2.0
	A_info.row[2]  = 2; A_info.col[2]  = 0; A_info.val[2]  = 2.0;
	// A_1,1  = 5.0
	A_info.row[3]  = 1; A_info.col[3]  = 1; A_info.val[3]  = 5.0;
	// A_2,1  = 1.0
	A_info.row[4]  = 2; A_info.col[4]  = 1; A_info.val[4]  = 1.0;
	// A_2,2  = 5.0
	A_info.row[5]  = 2; A_info.col[5]  = 2; A_info.val[5]  = 5.0;
	//
	// A_3,3  = 5.0
	A_info.row[6]  = 3; A_info.col[6]  = 3; A_info.val[6]  = 5.0;
	// A_4,3  = 4.0
	A_info.row[7]  = 4; A_info.col[7]  = 3; A_info.val[7]  = 4.0;
	// A_5,3  = 2.0
	A_info.row[8]  = 5; A_info.col[8]  = 3; A_info.val[8]  = 2.0;
	// A_4,4  = 5.0
	A_info.row[9]  = 4; A_info.col[9]  = 4; A_info.val[9]  = 5.0;
	// A_5,4 = 1.0
	A_info.row[10] = 5; A_info.col[10] = 4; A_info.val[10] = 1.0;
	// A_5,5 = 5.0
	A_info.row[11] = 5; A_info.col[11] = 5; A_info.val[11] = 5.0;

	// initialize the matrix using only the sparsity pattern
	cholmod_obj.init(A_info);

	// factor the matrix using the values
	cholmod_obj.update(A_info);

	// compute log of determinant of A
	double logdet_A = cholmod_obj.logdet();

	// check its value
	ok &= std::fabs( logdet_A / (2.0 * std::log(36.0)) - 1.0 ) <= eps;

	// test solve
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
		cholmod_obj.solve(row, val_in, val_out);
		//
		for(size_t k = 0; k < row.size(); k++)
		{	size_t i       = row[k];
			double check_i = A_inv[ i * nrow + j ];
			ok &= std::fabs( val_out[k] / check_i - 1.0 ) <= eps;
		}
	}
	return ok;
}
// END C++
