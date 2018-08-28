# ifndef CPPAD_MIXED_SPARSE_RCV2EIGEN_HPP
# define CPPAD_MIXED_SPARSE_RCV2EIGEN_HPP
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin sparse_rcv2eigen$$
$spell
	Eigen
	CppAD
	rcv
$$

$section Convert a CppAD Sparse Matrix to an Eigen Sparse Matrix$$

$head Syntax$$
$icode%m_eigen% = CppAD::mixed::sparse_rcv2eigen(%m_rcv%)%$$

$head Prototype$$
$srcfile%include/cppad/mixed/sparse_rcv2eigen.hpp
	%0%// BEGIN_PROTOTYPE%// END_PROTOTYPE%1
%$$

$head Private$$
This routine is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Scalar$$
is the element type for the sparse matrix.

$head m_rcv$$
is a CppAD representation of the sparse matrix.

$head m_eigen$$
The return value is an eigen representation of the sparse matrix.

$children%example/private/sparse_rcv2eigen.cpp
%$$
$head Example$$
The file $cref sparse_rcv2eigen.cpp$$ is an example
and test of $code sparse_rcv2eigen$$.

$end
-----------------------------------------------------------------------------
*/

# include <Eigen/SparseCore>
# include <cppad/mixed/typedef.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN_PROTOTYPE
template <class Scalar>
Eigen::SparseMatrix<Scalar> sparse_rcv2eigen(
	const CppAD::sparse_rcv<s_vector, CppAD::vector<Scalar> >& m_rcv
)
// END_PROTOTYPE
{	// m_rcv information
	size_t nr  = m_rcv.nr();
	size_t nc  = m_rcv.nc();
	size_t nnz = m_rcv.nnz();
	const s_vector&              row( m_rcv.row() );
	const s_vector&              col( m_rcv.col() );
	const CppAD::vector<Scalar>& val( m_rcv.val() );
	//
	// m_triplet corresponding to m_rcv
	typedef Eigen::Triplet<Scalar> triplet;
	std::vector<triplet> m_triplet(nnz);
	for(size_t k = 0; k < nnz; ++k)
		m_triplet[k] = triplet( int(row[k]), int(col[k]), val[k] );
	//
	// m_eigen corresponding to m_triplet
	Eigen::SparseMatrix<Scalar> m_eigen;
	m_eigen.resize( int(nr), int(nc) );
	m_eigen.setFromTriplets( m_triplet.begin(), m_triplet.end() );
	//
	return m_eigen;
}

} } // END_CPPAD_MIXED_NAMESPACE

# endif
