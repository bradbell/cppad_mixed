// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_SPARSE_INFO2EIGEN_HPP
# define CPPAD_MIXED_SPARSE_INFO2EIGEN_HPP

/*
{xrst_begin sparse_info2eigen dev}
{xrst_spell
   nc
   nr
}

Convert a sparse_mat_info Representation to An Eigen Sparse Matrix
##################################################################

Syntax
******
``CppAD::mixed::sparse_info2eigen`` ( *matrix* , *info* , *nr* , *nc* )

Private
*******
This routine is an implementation detail and not part of the
CppAD Mixed user API.

matrix
******
This argument has prototype

   ``const Eigen::SparseMatrix<double`` , *Option* , *Index* >& *matrix*

The input value of *matrix* does not matter.
Upon return, it contains the *nr* by *nc* matrix specified by
*info* .

Option
======
This must be ``Eigen::ColMajor`` .

Index
=====
The index type fro this sparse matrix are arbitrary; i.e.,
has no restrictions.

info
****
This argument has prototype

   ``const CppAD::mixed::sparse_mat_info&`` *info*

This object must be in
:ref:`column major<sparse_mat_info@Notation@Column Major Order>` order.

*nr*
This argument has prototype

   *size_t* ``nr``

and is the number of rows in the matrix.

*nc*
This argument has prototype

   *size_t* ``nc``

and is the number of columns in the matrix.
{xrst_toc_hidden
   example/private/sparse_info2eigen.cpp
}
Example
*******
The file :ref:`sparse_info2eigen.cpp-name` is an example
and test of ``sparse_info2eigen`` .

{xrst_end sparse_info2eigen}
*/
# include <Eigen/SparseCore>
# include <cppad/mixed/sparse_mat_info.hpp>

namespace CppAD { namespace mixed {

      template <class Index>
      void sparse_info2eigen(
         Eigen::SparseMatrix<double, Eigen::ColMajor, Index>& matrix ,
         const sparse_mat_info                                info   ,
         size_t                                               nr     ,
         size_t                                               nc     )
      {  using Eigen::ColMajor;
         typedef typename
         Eigen::SparseMatrix<double, ColMajor, Index>::Index  index;
         //
         size_t total_nnz = info.col.size();
         assert( info.row.size() == total_nnz );
         //
         // compute number of non-zeros in each column
         Eigen::Matrix<index, Eigen::Dynamic, 1> nnz(nc);
         size_t k = 0;
         for(size_t j = 0; j < nc; j++)
         {  nnz[j] = 0;
            while( k < total_nnz && info.col[k] == j )
            {  nnz[j]++;
               k++;
            }
         }
         // make sure info is in column major format
         assert( k == total_nnz );
         //
         // reserve space for the matrix
         matrix.resize( int(nr), int(nc) );
         matrix.reserve(nnz);
         //
         // set the values in the matrix
         for(k = 0; k < total_nnz; k++)
         {  size_t r = info.row[k];
            size_t c = info.col[k];
            matrix.insert( int(r), int(c) ) = info.val[k];
         }
         //
         return;
      }
} }
# endif
