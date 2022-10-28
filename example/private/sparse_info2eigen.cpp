// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
$begin sparse_info2eigen.cpp$$
$spell
	tri
	cppad
	sym
	eigen
$$

$section sparse_info2eigen: Example and Test$$

$head Private$$
This example is not part of the
$cref/cppad_mixed public API/base_class/$$.

$code
$srcthisfile%0%// BEGIN C++%// END C++%1%$$
$$
$end
*/
// BEGIN C++
# include <cppad/mixed/sparse_info2eigen.hpp>

bool sparse_info2eigen_xam(void)
{	bool ok = true;
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor, int> sparse_matrix;
	typedef sparse_matrix::InnerIterator                      iterator;
	//
	size_t nr = 3;
	size_t nc = 4;
	CppAD::mixed::sparse_mat_info info;
	size_t count = 0;
	// info in column major order
	for(size_t j = 0; j < nc; j++)
	{	for(size_t i = j; i < nr; i++)
		{	info.row.push_back(i);
			info.col.push_back(j);
			info.val.push_back( double(++count) );
		}
	}
	sparse_matrix matrix;
	CppAD::mixed::sparse_info2eigen(matrix, info, nr, nc);
	//
	// check the result values
	ok   &= matrix.outerSize() == matrix.cols();
	ok   &= size_t( matrix.rows() ) == nr;
	ok   &= size_t( matrix.cols() ) == nc;
	count = 0;
	for(size_t j = 0; j < nc; j++)
	{	for(iterator itr(matrix, int(j)); itr; ++itr)
		{	ok &= size_t( itr.row() )   == info.row[count];
			ok &= size_t( itr.col() )   == info.col[count];
			ok &= itr.value()           == info.val[count];
			++count;
		}
	}
	ok &= count == (nr * (nr + 1)) / 2;
	//
	return ok;
}
// END C++
