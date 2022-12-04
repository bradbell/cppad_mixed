// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_SPARSE_MAT2LOW_HPP
# define CPPAD_MIXED_SPARSE_MAT2LOW_HPP

/*
{xrst_begin sparse_mat2low}
{xrst_spell
   cols
}

Extract the Lower Triangular From an Eigen Symmetric Matrix
###########################################################

Syntax
******
*lower* = *CppAD::mixed::sparse_mat2low* ( ``matrix`` )

Private
*******
This routine is an implementation detail and not part of the
CppAD Mixed user API.

matrix
******
The argument has prototype

   ``const Eigen::SparseMatrix<`` *Scalar* , *Options* , *Index* >& *matrix*

and has the same number of rows as columns; i.e.

   *matrix* . ``rows`` () == *matrix* . ``cols`` ()

Scalar
======
The scalar type fro this sparse matrix are arbitrary; i.e,
has no restrictions.

Options
=======
The options for this sparse matrix are arbitrary; i.e.,
has no restrictions.

Index
=====
The index type fro this sparse matrix are arbitrary; i.e.,
has no restrictions.

lower
*****
The return value has prototype

   ``Eigen::SparseMatrix<`` *Scalar* , *Options* , *Index* > *lower*

and is the lower triangle of *matrix* ; i.e.,
it has the same lower triangle as *matrix* and it has no entries
above the diagonal.
{xrst_toc_hidden
   example/private/sparse_mat2low.cpp
}
Example
*******
The file :ref:`sparse_mat2low.cpp-name` is an example
and test of ``sparse_mat2low`` .

{xrst_end sparse_mat2low}
*/
# include <Eigen/SparseCore>

namespace CppAD { namespace mixed {

   template <class Sparse_matrix>
   Sparse_matrix sparse_mat2low(const Sparse_matrix& matrix)
   {  typedef typename Sparse_matrix::Index         index;
      typedef typename Sparse_matrix::InnerIterator iterator;
      assert( matrix.rows() == matrix.cols() );
      assert( matrix.rows() == matrix.outerSize() );
      //
      // determine the number of non-zeros entries for each outer iteration
      Eigen::Matrix<index, Eigen::Dynamic, 1> nnz( matrix.outerSize() );
      for(index k = 0; k < matrix.outerSize(); ++k)
      {  nnz[k] = index(0);
         for(iterator itr(matrix, k); itr; ++itr)
         {  if( itr.col() <= itr.row() )
               ++nnz[k]; // entries is at or below the diagonal
         }
      }
      //
      // reserve space for the result
      Sparse_matrix result( matrix.rows(), matrix.rows() );
      result.reserve(nnz);
      //
      for(index k = 0; k < matrix.outerSize(); ++k)
      {  for(iterator itr(matrix, k); itr; ++itr)
         {  if( itr.col() <= itr.row() )
               result.insert(itr.row(), itr.col()) = itr.value();
         }
      }
      return result;
   }
} }

# endif
