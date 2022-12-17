/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin sparse_eigen2rcv.cpp$$
$spell
	cppad
	eigen
	rcv
$$

$nospell
$bold This is old cppad_mixed documentation:$$ Here is a link to its
$href%http://bradbell.github.io/cppad_mixed%current documentation%$$.
$$
$section sparse_eigen2rcv: Example and Test$$

$head Private$$
This example is not part of the
$cref/cppad_mixed public API/base_class/$$.

$code
$srcthisfile%0%// BEGIN C++%// END C++%1%$$
$$
$end
*/
// BEGIN C++
# include <cppad/mixed/sparse_eigen2rcv.hpp>

bool sparse_eigen2rcv_xam(void)
{	bool ok = true;
	size_t nr = 3;
	size_t nc = nr;
	// declaring m_eigen( int(nr), int(nc) ) does not work.
	Eigen::SparseMatrix<double> m_eigen;
	//
	// m_eigen, nnz
	size_t nnz = 0;
	m_eigen.resize( int(nr), int(nc));
	for(size_t i = 0; i < nr; i++)
	{	for(size_t j = i; j < nc; j++)
		{	// nnz is index of this entry in row major order
			m_eigen.insert( int(i), int(j) ) = double(nnz);
			++nnz;
		}
	}
	// m_rcv
	CppAD::mixed::d_sparse_rcv m_rcv;
	m_rcv = CppAD::mixed::sparse_eigen2rcv(m_eigen);
	//
	// check sizes
	ok &= m_rcv.nr()  == nr;
	ok &= m_rcv.nc()  == nc;
	ok &= m_rcv.nnz() == nnz;
	//
	// check the result values
	const CppAD::mixed::s_vector& row( m_rcv.row() );
	const CppAD::mixed::s_vector& col( m_rcv.col() );
	const CppAD::mixed::d_vector& val( m_rcv.val() );
	const CppAD::mixed::s_vector  row_major = m_rcv.row_major();
	for(size_t k = 0; k < nnz; ++k)
	{	size_t i = row[ row_major[k] ];
		size_t j = col[ row_major[k] ];
		ok &= i <= j;
		ok &= val[ row_major[k] ] == double(k);
	}
	return ok;
}
// END C++
