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
$begin ldlt_cholmod_xam.cpp$$
$spell
	nrow
	init
	CppAD
	cholmod ldlt_obj
	logdet
	sim_cov
$$

$section Example Using Cholmod LDLT Factorization$$

$head Problem Description$$
We define the lower triangular matrix
$latex \[
	L =
	\left( \begin{array}{ccc}
		1 & 0 & 0 \\
		1 & 2 & 0 \\
		1 & 2 & 3
	\end{array} \right)
	\W{,}
	L^\R{T} =
	\left( \begin{array}{ccc}
		1 & 1 & 1 \\
		0 & 2 & 2 \\
		0 & 0 & 3
	\end{array} \right)
\] $$
and the positive definite matrix
$latex \[
	H = L L^\R{T} =
	\left( \begin{array}{ccc}
		1 & 1 & 1 \\
		1 & 5 & 5 \\
		1 & 5 & 14
	\end{array} \right)
\] $$
The inverse of $latex H$$ is given by
$latex \[
	H^{-1} = L^\R{-T} L^{-1} =
	\frac{1}{36}
	\left( \begin{array}{ccc}
		45  & -9  & 0  \\
		-9  & 13  & -4 \\
		0   & -4  & 4
	\end{array} \right)
\] $$
which can be checked by multiplying by $latex H H^{-1}$$.

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
	logdet_H = ldlt_obj.logdet(negative);
$$

$head solve_H$$
See the following under Source Code below:
$codep
	ldlt_obj.solve_H(row, val_in, val_out);
$$

$head sim_cov$$
See the following under Source Code below:
$codep
	ok &= ldlt_obj.sim_cov(w, v)
$$

$head Source Code$$
$code
$srcfile%example/private/ldlt_cholmod_xam.cpp%5%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/mixed/ldlt_cholmod.hpp>
# include <limits>
# include <cmath>
# include <cassert>

bool ldlt_cholmod_xam(void)
{	bool ok    = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();

	double H_inv[] = {
		45.0,  -9.0,  0.0,
		-9.0,  13.0, -4.0,
		0.0,   -4.0,  4.0
	};
	for(size_t i = 0; i < sizeof(H_inv)/sizeof(H_inv[0]); i++)
		H_inv[i] /= 36.;

	// create cholmod object
	size_t nrow = 3;    // number of rows in H
	size_t ncol = nrow; // number of columns in H
	CppAD::mixed::ldlt_cholmod ldlt_obj(nrow);
	assert( nrow * ncol == sizeof(H_inv) / sizeof(H_inv[0]) );

	// create a sparse matrix representation of the lower triangular of H
	CppAD::mixed::sparse_mat_info H_info;
	H_info.resize(6);
	// H_0,0  = 1.0
	H_info.row[0]  = 0; H_info.col[0]  = 0; H_info.val[0]  = 1.0;
	// H_1,0  = 1.0
	H_info.row[1]  = 1; H_info.col[1]  = 0; H_info.val[1]  = 1.0;
	// H_2,0  = 1.0
	H_info.row[2]  = 2; H_info.col[2]  = 0; H_info.val[2]  = 1.0;
	// H_1,1  = 5.0
	H_info.row[3]  = 1; H_info.col[3]  = 1; H_info.val[3]  = 5.0;
	// H_2,1  = 5.0
	H_info.row[4]  = 2; H_info.col[4]  = 1; H_info.val[4]  = 5.0;
	// H_2,2  = 14.0
	H_info.row[5]  = 2; H_info.col[5]  = 2; H_info.val[5]  = 14.0;
	//
	// initialize the matrix using only the sparsity pattern
	ldlt_obj.init(H_info);

	// factor the matrix using the values
	ldlt_obj.update(H_info);

	// compute log of determinant of H
	size_t negative;
	double logdet_H = ldlt_obj.logdet(negative);
	ok &= negative == 0;

	// check its value
	ok &= std::fabs( logdet_H / std::log(36.0) - 1.0 ) <= eps;

	// test solve_H
	{	CppAD::vector<size_t> row(3);
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
	}

	// test inv
	{	size_t K = 9;
		CppAD::vector<size_t> row_in(K), col_in(K);
		CppAD::vector<double> val_out(K);
		{	size_t k = 0;
			// use row major ordering to test sorting
			for(size_t i = 0; i < nrow; i++)
			{	for(size_t j = 0; j < nrow; j++)
				{	row_in[k] = i;
					col_in[k] = j;
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
	//
	// test sim_cov
	CppAD::vector<double> w(3), v(3), c(3);
	for(size_t i = 0; i < 3; i++)
		w[i] = double( 2 * i + 1);
	ok &= ldlt_obj.sim_cov(w, v);
	// solve w = L^{T} c
	c[2] = w[2] / 3.0;
	c[1] = ( w[1] - 2 * c[2] ) / 2.0;
	c[0] = ( w[0] - 1.0 * c[1] - 1.0 * c[2] ) / 1.0;
	for(size_t i = 0; i < 3; i++)
		ok  &= std::fabs( v[i] - c[i] ) <= eps;

	return ok;
}
// END C++
