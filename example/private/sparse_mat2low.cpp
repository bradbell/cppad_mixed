// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-20 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin sparse_mat2low.cpp$$
$spell
	tri
	cppad
	sym
$$

$section sparse_mat2low: Example and Test$$

$head Private$$
This example is not part of the
$cref/cppad_mixed public API/base_class/$$.

$code
$srcthisfile%0%// BEGIN C++%// END C++%1%$$
$$
$end
*/
// BEGIN C++
# include <cstddef>
# include <cppad/mixed/sparse_mat2low.hpp>

bool sparse_mat2low_xam(void)
{	bool ok = true;
	typedef Eigen::SparseMatrix<size_t> sparse_matrix;
	//
	size_t nr = 3;
	sparse_matrix matrix;
	matrix.resize(int(nr), int(nr));
	Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> check(nr, nr);
	size_t count = 0;
	for(size_t i = 0; i < nr; i++)
	{	for(size_t j = 0; j < nr; j++)
		{	matrix.insert(int(i), int(j)) = count;
			check(i, j)                   = count;
			++count;
		}
	}
	sparse_matrix lower = CppAD::mixed::sparse_mat2low(matrix);
	//
	// check the result values
	count = 0;
	for(size_t k = 0; k < nr; k++)
	{	for(sparse_matrix::InnerIterator itr(lower, int(k)); itr; ++itr)
		{	int i = itr.row();
			int j = itr.col();
			ok   &= itr.value() == check(i, j);
			// make sure only lower triangle is inclued
			ok   &= j <= i;
			++count;
		}
	}
	// make sure all of lower traingle is included
	ok &= count == (nr * (nr + 1) ) / 2;
	return ok;
}
// END C++
