// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-23 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_SPARSE_SCALE_DIAG_HPP
# define CPPAD_MIXED_SPARSE_SCALE_DIAG_HPP

/*
{xrst_begin sparse_scale_diag dev}
{xrst_spell
   diag
}

Scales the Diagonal of an Eigen Sparse Matrix
#############################################

Syntax
******
``CppAD::mixed::sparse_scale_diag`` ( *scale* , *matrix* )

Purpose
*******
Each of the diagonal elements of *matrix* is replaced by
*scale* times its original value.

Private
*******
This routine is an implementation detail and not part of the
CppAD Mixed user API.

scale
*****
This argument has prototype

   ``const`` *Float* & *scale*

where the type *Float* can be converted to the scalar
type for the matrix.

matrix
******
This argument has prototype

   ``Eigen::SparseMatrix<`` *Scalar* , *Options* , *Index* >& *matrix*

This is the sparse matrix for which we are scaling the diagonal elements.
The sparsity pattern for the matrix is not modified.

Scalar
======
This type must support the conversion

   ``static_cast`` < *Scalar* >( *scale* )

.

Options
=======
The options for this sparse matrix are arbitrary; i.e.,
has no restrictions.

Index
=====
The index type fro this sparse matrix are arbitrary; i.e.,
has no restrictions.
{xrst_toc_hidden
   example/private/sparse_scale_diag.cpp
}
Example
*******
The file :ref:`sparse_scale_diag.cpp-name` is an example
and test of ``sparse_scale_diag`` .

{xrst_end sparse_scale_diag}
*/
# include <Eigen/SparseCore>

namespace CppAD { namespace mixed {
   template <class Float, class Sparse_matrix>
   void sparse_scale_diag(Float scale, Sparse_matrix& matrix)
   {  typedef typename Sparse_matrix::Index         index;
      typedef typename Sparse_matrix::Scalar        scalar;
      typedef typename Sparse_matrix::InnerIterator iterator;
      scalar factor = static_cast<scalar>(scale);
      for(index k = 0; k < matrix.outerSize(); ++k)
      {  for(iterator itr(matrix, k); itr; ++itr)
         {  if( itr.row() == itr.col() )
               itr.valueRef() = itr.value() * factor;
         }
      }
   }
} }

# endif
