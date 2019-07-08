# ifndef CPPAD_MIXED_SPARSE_EIGEN2RCV_HPP
# define CPPAD_MIXED_SPARSE_EIGEN2RCV_HPP
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-19 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin sparse_eigen2rcv$$
$spell
	Eigen
	CppAD
	rcv
	Cppad
$$

$section Convert an Eigen Sparse Matrix to a Cppad Sparse Matrix$$

$head Syntax$$
$icode%m_rcv% = CppAD::mixed::sparse_eigen2rcv(%m_eigen%)%$$

$head Prototype$$
$srcfile%include/cppad/mixed/sparse_eigen2rcv.hpp
	%0%// BEGIN_PROTOTYPE%// END_PROTOTYPE%1
%$$

$head Private$$
This routine is an implementation detail and not part of the
CppAD Mixed user API.

$head Scalar$$
is the element type for the sparse matrix.

$head m_eigen$$
is an Eigen representation of the sparse matrix.

$head m_rcv$$
is a Cppad representation of the matrix.

$children%example/private/sparse_eigen2rcv.cpp
%$$
$head Example$$
The file $cref sparse_eigen2rcv.cpp$$ is an example
and test of $code sparse_eigen2rcv$$.

$end
-----------------------------------------------------------------------------
*/

# include <Eigen/SparseCore>
# include <cppad/mixed/typedef.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN_PROTOTYPE
template <class Scalar>
CppAD::sparse_rcv<s_vector, CppAD::vector<Scalar> > sparse_eigen2rcv(
	const Eigen::SparseMatrix<Scalar>& m_eigen
)
// END_PROTOTYPE
{	//
	// iterator
	typedef typename Eigen::SparseMatrix<Scalar>::InnerIterator iterator;
	//
	// row, col, val, nnz
	s_vector              row;
	s_vector              col;
	CppAD::vector<Scalar> val;
	size_t nnz = 0;
	for(int j = 0; j < m_eigen.outerSize(); ++j)
	{	for(iterator itr(m_eigen, j); itr; ++itr)
		{	row.push_back( size_t(itr.row()) );
			col.push_back( size_t(itr.col()) );
			val.push_back( itr.value() );
			++nnz;
		}
	}
	//
	// m_rc
	size_t    nr  = size_t( m_eigen.rows() );
	size_t    nc  = size_t( m_eigen.cols() );
	sparse_rc m_rc(nr, nc, nnz);
	for(size_t k = 0; k < nnz; ++k)
		m_rc.set(k, row[k], col[k]);
	//
	// m_rcv
	CppAD::sparse_rcv<s_vector, CppAD::vector<Scalar> > m_rcv(m_rc);
	for(size_t k = 0; k < nnz; ++k)
		m_rcv.set(k, val[k]);
	//
	return m_rcv;
}

} } // END_CPPAD_MIXED_NAMESPACE

# endif
