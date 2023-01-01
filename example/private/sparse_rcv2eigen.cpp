/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin sparse_rcv2eigen.cpp$$
$spell
	cppad
	eigen
	rcv
$$

$nospell
$bold This is cppad_mixed--20220519 documentation:$$ Here is a link to its
$href%https://cppad-mixed.readthedocs.io%current documentation%$$.
$$
$section sparse_rcv2eigen: Example and Test$$

$head Private$$
This example is not part of the
$cref/cppad_mixed public API/base_class/$$.

$code
$srcthisfile%0%// BEGIN C++%// END C++%1%$$
$$
$end
*/
// BEGIN C++
# include <cppad/mixed/sparse_rcv2eigen.hpp>

bool sparse_rcv2eigen_xam(void)
{	bool ok = true;
	size_t nr = 3;
	size_t nc = nr;
	size_t nnz = (nr * (nr + 1)) / 2;
	//
	// m_rc
	CppAD::mixed::sparse_rc m_rc(nr, nc, nnz);
	size_t k = 0;
	for(size_t i = 0; i < nr; i++)
	{	for(size_t j = i; j < nc; j++)
			m_rc.set(k++, i, j);
	}
	ok &= k == nnz;
	//
	// m_rcv
	CppAD::mixed::d_sparse_rcv m_rcv( m_rc );
	for(k = 0; k < nnz; ++k)
		m_rcv.set(k, double(k));
	//
	Eigen::SparseMatrix<double, Eigen::RowMajor> m_eigen;
	m_eigen = CppAD::mixed::sparse_rcv2eigen(m_rcv);
	//
	// check sizes
	ok &= m_eigen.rows() == int(nr);
	ok &= m_eigen.cols() == int(nc);
	//
	// check the result values
	typedef typename
	Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator iterator;
	k = 0;
	for(int i = 0; i < int(nr); ++i)
	{	for(iterator itr(m_eigen, i); itr; ++itr)
		{	ok &= itr.value() == double(k++);
			ok &= itr.row() == i;
			ok &= i <= itr.col();
		}
	}
	ok &= k == nnz;
	//
	return ok;
}
// END C++
