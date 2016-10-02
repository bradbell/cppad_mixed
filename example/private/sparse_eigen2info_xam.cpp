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
$begin sparse_eigen2info_xam.cpp$$
$spell
	tri
	cppad
	sym
	eigen
$$

$section sparse_eigen2info: Example and Test$$

$head Private$$
This example is not part of the
$cref/cppad_mixed public API/public/$$.

$code
$srcfile%example/private/sparse_eigen2info_xam.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$
$end
*/
// BEGIN C++
# include <cppad/mixed/sparse_eigen2info.hpp>

bool sparse_eigen2info_xam(void)
{	bool ok = true;
	typedef Eigen::SparseMatrix<double, Eigen::RowMajor, int> sparse_matrix;
	//
	int nr = 3;
	sparse_matrix matrix(nr, nr);
	size_t count = 0;
	for(int i = 0; i < nr; i++)
	{	for(int j = 0; j <= i; j++)
		{	matrix.insert(i, j) = double(count);
			++count;
		}
	}
	CppAD::mixed::sparse_mat_info info =
		CppAD::mixed::sparse_eigen2info(matrix);
	//
	// check the result values
	count = 0;
	for(size_t i = 0; i < size_t(nr); i++)
	{	for(size_t j = 0; j <= i; j++)
		{	ok &= info.row[count] == i;
			ok &= info.col[count] == j;
			ok &= info.val[count] == double(count);
			++count;
		}
	}
	ok &= count == info.row.size();
	ok &= count == info.col.size();
	ok &= count == info.val.size();
	//
	return ok;
}
// END C++
