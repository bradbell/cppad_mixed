// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_SPARSE_LOW2SYM_HPP
# define CPPAD_MIXED_SPARSE_LOW2SYM_HPP

/*
{xrst_begin sparse_low2sym}
{xrst_spell
   cols
   sym
}

Convert an Eigen Lower Triangular Matrix To a Symmetric Matrix
##############################################################

Syntax
******
*symmetric* = *CppAD::mixed::sparse_low2sym* ( ``lower`` )

Private
*******
This routine is an implementation detail and not part of the
CppAD Mixed user API.

lower
*****
The argument has prototype

   ``const Eigen::SparseMatrix<`` *Scalar* , *Options* , *Index* >& *lower*

and has the same number of rows as columns; i.e.

   ``lower`` . *rows* () == ``lower.cols`` ()

This is a lower triangular sparse matrix; i.e.,
for each entry in *lower* has a row index
that is greater than or equal its column index.

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

symmetric
*********
The return value has prototype

   ``Eigen::SparseMatrix<`` *Scalar* , *Options* , *Index* > *symmetric*

Its lower triangle has the same entries as *lower*
and it is a symmetric matrix; i.e.,
the entries above the diagonal have been set using the corresponding
entry below the diagonal.
{xrst_toc_hidden
   example/private/sparse_low2sym.cpp
}
Example
*******
The file :ref:`sparse_low2sym.cpp-name` is an example
and test of ``sparse_low2sym`` .

{xrst_end sparse_low2sym}
*/
# include <Eigen/SparseCore>

namespace CppAD { namespace mixed {
   template <class Sparse_matrix>
   Sparse_matrix sparse_low2sym(const Sparse_matrix& lower)
   {  typedef typename Sparse_matrix::Index         index;
      typedef typename Sparse_matrix::InnerIterator iterator;
      assert( lower.rows() == lower.cols() );
      assert( lower.rows() == lower.outerSize() );
      //
      // determine the number of non-zeros in each row or column of result
      Eigen::Matrix<index, Eigen::Dynamic, 1> nnz( lower.outerSize() );
      for(index k = 0; k < lower.outerSize(); ++k)
      {  nnz[k] = index(0);
         for(iterator itr(lower, k); itr; ++itr)
         {  ++nnz[k]; // cound entries in lower triangle
            if( itr.row() != itr.col() )
               ++nnz[k]; // entries above the diagonal
         }
      }
      //
      // reserve space for the result
      Sparse_matrix result( lower.rows(), lower.rows() );
      result.reserve(nnz);
      //
      for(index k = 0; k < lower.outerSize(); ++k)
      {  for(iterator itr(lower, k); itr; ++itr)
         {  index i = itr.row();
            index j = itr.col();
            assert( j <= i );
            result.insert(i, j) = itr.value();
            if( i != j )
               result.insert(j, i) = itr.value();
         }
      }
      return result;
   }
} }

# endif
