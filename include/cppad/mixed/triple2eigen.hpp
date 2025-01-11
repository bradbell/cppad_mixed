# ifndef CPPAD_MIXED_TRIPLE2EIGEN_HPP
# define CPPAD_MIXED_TRIPLE2EIGEN_HPP
// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin triple2eigen dev}
{xrst_spell
  nc
  nr
}

Convert Row, Column, Value Triple to an Eigen Sparse Matrix
###########################################################

Syntax
******

| ``CppAD::mixed::triple2eigen`` (
| |tab| *mat* , *nr* , *nc* , *row* , *col* , *val*
| )

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
If *d* is a ``double`` value,  *Scalar* ( ``d`` )
must be a corresponding element of the matrix.

nr
**
is the number of rows in the matrix.

nc
**
is the number of columns in the matrix.

row
***
contains the row indices for possibly non-zero elements of the matrix.

col
***
has the same size as *row* and
contains the column indices for possibly non-zero elements of the matrix.

val
***

Sparsity Pattern
================
If *val* . ``size`` () == 0 ,
the values in the matrix are not specified.
To be specific, for *k* = 0 , ..., *row* . ``size`` () ``-1`` ,
the element with index
( *row* [ *k* ], *col* [ *k* ]) has an unspecified value.

Sparse Matrix
=============
If *val* . ``size`` () != 0 ,
it contains the possibly non-zero values in the matrix.
To be specific, for *k* = 0 , ..., *row* . ``size`` () ``-1`` ,
the element with index ( *row* [ *k* ], *col* [ *k* ]) has value
*val* [ *k* ] .

mat
***
is a sparse representation of the specified matrix.
The input size and values in the matrix do not matter.
Upon return it is a sparse matrix with the specified size
and element values.

{xrst_end triple2eigen}
-----------------------------------------------------------------------------
*/

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

template <class Scalar>
// BEGIN_PROTOTYPE
void triple2eigen(
   Eigen::SparseMatrix<Scalar>&  mat  ,
   size_t                        nr   ,
   size_t                        nc   ,
   const s_vector&               row  ,
   const s_vector&               col  ,
   const CppAD::vector<Scalar>&  val  )
// END_PROTOTYPE
{  assert( row.size() == col.size() );
   assert( val.size() ==  0 || row.size() == val.size() );
   //
   mat.resize( int(nr), int(nc) );
   if( val.size() == 0 )
   {  for(size_t k = 0; k < row.size(); k++)
         mat.insert( int(row[k]), int(col[k]) ) = Scalar(0.0);
   }
   else
   {  for(size_t k = 0; k < row.size(); k++)
         mat.insert( int(row[k]), int(col[k]) ) = val[k];
   }
   return;
}

} } // END_CPPAD_MIXED_NAMESPACE

# endif
