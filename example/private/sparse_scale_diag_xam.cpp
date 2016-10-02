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
$begin sparse_scale_diag_xam.cpp$$
$spell
	tri
	cppad
$$

$section sparse_scale_diag: Example and Test$$

$head Private$$
This example is not part of the
$cref/cppad_mixed public API/public/$$.

$head Description$$
The call to $cref sparse_scale_diag$$ below computes the matrix
$codei%
	[ 3*1 0   0   ]
	[ 2   3*3 0   ]
	[ 4   5   3*6 ]
%$$

$code
$srcfile%example/private/sparse_scale_diag_xam.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$
$end
*/
// BEGIN C++
# include <cmath>
# include <cstddef>
# include <limits>
# include <cppad/mixed/sparse_scale_diag.hpp>

bool sparse_scale_diag_xam(void)
{	bool ok = true;
	double eps = 10. * std::numeric_limits<double>::epsilon();
	typedef Eigen::SparseMatrix<double> sparse_matrix;
	//
	int nr = 3;
	sparse_matrix matrix(nr, nr);
	int count = 0;
	for(int i = 0; i < nr; i++)
	{	for(int j = 0; j <= i; j++)
			matrix.insert(i, j) = double( ++count );
	}
	int scale = 3;
	CppAD::mixed::sparse_scale_diag(scale, matrix);
	//
	// Check that the diagonal has been scaled by there and that the
	// matrx is lower triangular
	for(int k = 0; k < nr; k++)
	{	for(sparse_matrix::InnerIterator itr(matrix, k); itr; ++itr)
		{	int i = itr.row();
			int j = itr.col();
			ok   &= j <= i;
			if( i == j )
			{	int original_value = (i + 1) * (i + 2) / 2;
				int new_value      = scale * original_value;
				ok &= std::fabs( itr.value() - new_value ) < eps;
			}
		}
	}
	return ok;
}
// END C++
