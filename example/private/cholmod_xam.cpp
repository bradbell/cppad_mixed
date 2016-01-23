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
We are given the matrix
$latex \[
	A = \left( \begin{array}{ccc}
		5 & 4 & 2 \\
		4 & 5 & 1 \\
		2 & 1 & 5
	\end{array} \right)
\] $$
We use $latex A^k$$ to denote the upper-left $latex k \times k$$
principal minor. The determinant of its principal minors are:
$latex \[
\begin{array}{rcl}
	\det \left( A^1 \right) & = & 5 \\
	\det \left( A^2 \right) & = & 9 \\
	\det \left( A^3 \right) & = & 36
\end{array}
\] $$
In addition
$latex \[
	A^{-1} = \frac{1}{36}
	\left( \begin{array}{ccc}
		24  & -18 & -6 \\
		-18 & 21  & 3 \\
		-6  & 3   & 9
	\end{array} \right)
\] $$
which can be checked by multiplying by $latex A$$.


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
	cholmod_obj.solve(row_in, val_in, row_out, val_out);
$$

$head Source Code$$
$code
$verbatim%example/private/cholmod_xam.cpp%5%// BEGIN C++%// END C++%1%$$
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
		 24.0, -18.0, -6.0,
		-18.0,  21.0,  3.0,
		 -6.0,   3.0,  9.0
	};
	for(size_t i = 0; i < sizeof(A_inv)/sizeof(A_inv[0]); i++)
		A_inv[i] /= 36.;

	// create cholmod object
	size_t nrow = 3;   // number of rows in A
	size_t ncol = nrow;
	CppAD::mixed::cholmod cholmod_obj(nrow);

	// create a sparse matrix representation of the lower triangular of A
	CppAD::mixed::sparse_mat_info A_info;
	A_info.resize(6);
	// A_0,0 = 5.0
	A_info.row[0] = 0; A_info.col[0] = 0; A_info.val[0] = 5.0;
	// A_1,0 = 4.0
	A_info.row[1] = 1; A_info.col[1] = 0; A_info.val[1] = 4.0;
	// A_2,0 = 2.0
	A_info.row[2] = 2; A_info.col[2] = 0; A_info.val[2] = 2.0;
	// A_1,1 = 5.0
	A_info.row[3] = 1; A_info.col[3] = 1; A_info.val[3] = 5.0;
	// A_2,1 = 1.0
	A_info.row[4] = 2; A_info.col[4] = 1; A_info.val[4] = 1.0;
	// A_2,2 = 5.0
	A_info.row[5] = 2; A_info.col[5] = 2; A_info.val[5] = 5.0;

	// initialize the matrix using only the sparsity pattern
	cholmod_obj.init(A_info);

	// factor the matrix using the values
	cholmod_obj.update(A_info);

	// compute log of determinant of A
	double logdet_A = cholmod_obj.logdet();

	// check its value
	ok &= std::fabs( logdet_A / std::log(36.0) - 1.0 ) <= eps;

	CppAD::vector<size_t> row_in(1), row_out(nrow);
	CppAD::vector<double> val_in(1), val_out(nrow);
	for(size_t i = 0; i < nrow; i++)
		row_out[i] = i;
	for(size_t j = 0; j < ncol; j++)
	{	// solve for the j-th column of the inverse matrix
		row_in[0] = j;
		val_in[0] = 1.0;
		cholmod_obj.solve(row_in, val_in, row_out, val_out);
		for(size_t i = 0; i < nrow; i++)
			ok &= std::fabs( val_out[i] / A_inv[i*ncol+j] - 1.0 ) <= eps;
	}
	return ok;
}
// END C++
