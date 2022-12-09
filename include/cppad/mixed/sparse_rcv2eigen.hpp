# ifndef CPPAD_MIXED_SPARSE_RCV2EIGEN_HPP
# define CPPAD_MIXED_SPARSE_RCV2EIGEN_HPP
// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin sparse_rcv2eigen dev}

Convert a CppAD Sparse Matrix to an Eigen Sparse Matrix
#######################################################

Syntax
******
*m_eigen* = ``CppAD::mixed::sparse_rcv2eigen`` ( *m_rcv* )

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE
   // END_PROTOTYPE
}

Private
*******
This routine is an implementation detail and not part of the
CppAD Mixed user API.

Scalar
******
is the element type for the sparse matrix.

m_rcv
*****
is a CppAD representation of the sparse matrix.

m_eigen
*******
The return value is an eigen representation of the sparse matrix.
{xrst_toc_hidden
   example/private/sparse_rcv2eigen.cpp
}
Example
*******
The file :ref:`sparse_rcv2eigen.cpp-name` is an example
and test of ``sparse_rcv2eigen`` .

{xrst_end sparse_rcv2eigen}
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
{  // m_rcv information
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
