// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_SPARSE_EIGEN2INFO_HPP
# define CPPAD_MIXED_SPARSE_EIGEN2INFO_HPP

/*
{xrst_begin sparse_eigen2info}

Convert An Eigen Sparse Matrix to a sparse_mat_info Representation
##################################################################

Syntax
******
``CppAD::mixed::sparse_eigen2info`` ( *matrix* , *info* )

Private
*******
This routine is an implementation detail and not part of the
CppAD Mixed user API.

matrix
******
The argument has prototype

   ``const Eigen::SparseMatrix<double`` , *Option* , *Index* >& *matrix*

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

   ``CppAD::mixed::sparse_mat_info`` *info*

There are different cases depending on if it is the
:ref:`sparse_mat_info@Notation@Empty Matrix` on input;

Empty on Input
==============
Upon return *info* has the same sparse matrix information
as *matrix* and is in
:ref:`column major<sparse_mat_info@Notation@Column Major Order>` order.
In this case it is assumed that there is at least one entry in every
column of *matrix* .

Non-Empty on Input
==================
It is assumed that on input, the size of *info* . ``val``
is the same as *info* . ``row`` and *info* . ``col`` .
Upon return, the following conditions hold:

#. The sparsity pattern in *info* is not modified.
#. The value of any elements in *matrix* ,
   that also appear in the *info* , are copied to *info* .
#. Any elements that appear in *info* , but not in *matrix* ,
   have value zero in the return value of *info* .
#. The elements of *matrix* ,
   that do not appear in *info* , are ignored.

{xrst_toc_hidden
   example/private/sparse_eigen2info.cpp
}
Example
*******
The file :ref:`sparse_eigen2info.cpp-name` is an example
and test of ``sparse_eigen2info`` .

{xrst_end sparse_eigen2info}
*/
# include <Eigen/SparseCore>
# include <cppad/mixed/sparse_mat_info.hpp>

namespace CppAD { namespace mixed {

      template <class Index>
      void sparse_eigen2info(
         const Eigen::SparseMatrix<double, Eigen::ColMajor, Index>& matrix ,
         sparse_mat_info&                                           info   )
      {  using Eigen::ColMajor;
         typedef typename Eigen::
         SparseMatrix<double, ColMajor, Index>::InnerIterator iterator;
         typedef typename Eigen::
         SparseMatrix<double, ColMajor, Index>::Index         index;
         //
         // case where input value of info is empty matrix
         if( info.row.size() == 0 )
         {  for(index j = 0; j < matrix.outerSize(); j++)
            {  for(iterator itr(matrix, j); itr; ++itr)
               {  info.row.push_back( size_t( itr.row() ) );
                  info.col.push_back( size_t( itr.col() ) );
                  info.val.push_back( itr.value() );
               }
            }
            return;
         }
         // case where input value of info is non-empty
         //
         // initilize all the values as zero
         size_t K = info.row.size();
         for(size_t k = 0; k < K; k++)
            info.val[k] = 0.0;
         //
         size_t k = 0;
         for(index j = 0; j < matrix.outerSize(); j++)
         {  // skip entries in info that are not in matrix
            size_t c = size_t(j);
            while( info.col[k] < c )
               k++;
            assert( info.col[k] == c );
            for(iterator itr(matrix, j); itr; ++itr)
            {  // skip entries in info that are not in matrix
               size_t r = size_t( itr.row() );
               while( k < K && info.col[k] == c && info.row[k] < r )
                  k++;
               // check if this entry is in both matrix and info
               if( k < K && info.col[k] == c && info.row[k] == r )
                  info.val[k++] = itr.value();
            }
         }
         return;
      }
} }
# endif
