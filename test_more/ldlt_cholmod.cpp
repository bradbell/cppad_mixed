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
begin ldlt_cholmod_xam.cpp
$spell
	nrow
	init
	Cholesky
	CppAD
	cholmod ldlt_obj
	logdet
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
$cref/Source Code/ldlt_cholmod_xam.cpp/Source Code/$$ below:
$codep
	ldlt_obj.init(H_info);
$$

$head update$$
See the following under Source Code below:
$codep
	ldlt_obj.update(H_info);
$$

$head logdet$$
See the following under Source Code below:
$codep
	ldlt_obj.logdet(H_info);
$$

$head solve_H$$
See the following under Source Code below:
$codep
	ldlt_obj.solve_H(row, val_in, val_out);
$$

$head Source Code$$
$code
$srcfile%test_more/ldlt_cholmod.cpp%5%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/mixed/ldlt_cholmod.hpp>
# include <limits>
# include <cmath>
# include <cassert>

bool ldlt_cholmod(void)
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

	// create a sparse matrix representation of the lower triangular of H
	CppAD::mixed::sparse_mat_info H_info;
	H_info.resize(12);
	// H_0,0  = 5.0
	H_info.row[0]  = 0; H_info.col[0]  = 0; H_info.val[0]  = 5.0;
	// H_1,0  = 4.0
	H_info.row[1]  = 1; H_info.col[1]  = 0; H_info.val[1]  = 4.0;
	// H_2,0  = 2.0
	H_info.row[2]  = 2; H_info.col[2]  = 0; H_info.val[2]  = 2.0;
	// H_1,1  = 5.0
	H_info.row[3]  = 1; H_info.col[3]  = 1; H_info.val[3]  = 5.0;
	// H_2,1  = 1.0
	H_info.row[4]  = 2; H_info.col[4]  = 1; H_info.val[4]  = 1.0;
	// H_2,2  = 5.0
	H_info.row[5]  = 2; H_info.col[5]  = 2; H_info.val[5]  = 5.0;
	//
	// H_3,3  = 5.0
	H_info.row[6]  = 3; H_info.col[6]  = 3; H_info.val[6]  = 5.0;
	// H_4,3  = 4.0
	H_info.row[7]  = 4; H_info.col[7]  = 3; H_info.val[7]  = 4.0;
	// H_5,3  = 2.0
	H_info.row[8]  = 5; H_info.col[8]  = 3; H_info.val[8]  = 2.0;
	// H_4,4  = 5.0
	H_info.row[9]  = 4; H_info.col[9]  = 4; H_info.val[9]  = 5.0;
	// H_5,4 = 1.0
	H_info.row[10] = 5; H_info.col[10] = 4; H_info.val[10] = 1.0;
	// H_5,5 = 5.0
	H_info.row[11] = 5; H_info.col[11] = 5; H_info.val[11] = 5.0;

	// initialize the matrix using only the sparsity pattern
	ldlt_obj.init(H_info);

	// factor the matrix using the values
	ldlt_obj.update(H_info);

	// compute log of determinant of H
	int sign;
	double logdet_H = ldlt_obj.logdet(sign);
	ok &= sign == 1;

	// check its value
	ok &= std::fabs( logdet_H / (2.0 * std::log(36.0)) - 1.0 ) <= eps;

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
		ldlt_obj.solve_H(row, val_in, val_out);
		//
		for(size_t k = 0; k < row.size(); k++)
		{	size_t i       = row[k];
			double check_i = H_inv[ i * nrow + j ];
			ok &= std::fabs( val_out[k] / check_i - 1.0 ) <= eps;
		}
	}
	return ok;
}
// END C++
