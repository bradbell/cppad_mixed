// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_SPARSE_MAT_INFO_HPP
# define CPPAD_MIXED_SPARSE_MAT_INFO_HPP
/*
{xrst_begin sparse_mat_info}
{xrst_spell
  resize
}

Sparse Matrix Information
#########################

Syntax
******

   *CppAD::mixed::sparse_mat_info* ``mat_info``

*mat_info* . ``resize`` ( *size* )

Purpose
*******
This structure holds information about a sparse matrix.

row
***
The field *mat_info* . ``row`` has prototype

   ``CppAD::vector<size_t>`` *mat_info* . ``row``

It has size zero when it is constructed.
After initialization it should contain the row indices
corresponding to possibly non-zero elements of the matrix.

K
=
We use *K* = *mat_info* . ``row.size`` () below.

col
***
The field *mat_info* . ``col`` has prototype

   ``CppAD::vector<size_t>`` *mat_info* . ``col``

It has size zero when it is constructed.
After initialization it should have the same size as *row*
and contain the column indices
corresponding to possibly non-zero elements of the matrix.

val
***
The field *mat_info* . ``val`` has prototype

   ``CppAD::vector<double>`` *mat_info* . ``val``

It has size zero when it is constructed.
After initialization it should either have size zero,
or the same size as *row* .

resize
******
The ``resize`` argument has prototype

   ``size_t`` *size*

All of the vectors,
*row* , *col* , and *val* ,
are modified to have the specified size.

Notation
********

Sparsity Pattern
================
We say that *mat_info* is a sparsity pattern if,
for *k* = 0 , ... , *K* ``-1`` ,
the element with index

   ( *mat_info* . ``row`` [ *k* ], *mat_info* . ``col`` [ *k* ])

is possibly non-zero and the size or elements of
*mat_info* . ``val`` are not specified.

Sparse Matrix
=============
We say that *mat_info* is a sparse matrix if,
for *k* = 0 , ... , *K* ``-1`` ,
the element with index

   ( *mat_info* . ``row`` [ *k* ], *mat_info* . ``col`` [ *k* ])

is possibly non-zero and has value *mat_info* . ``val`` [ *k* ] .

Empty Matrix
============
If *K* is zero ( *mat_info* . ``row.size`` () is zero),
we say that *mat_info* is the empty matrix.

Column Major Order
==================
If for *k* = 0 , ... , *K* ``-1`` ,

| |tab| *mat_info* . ``col`` [ *k* ] <= *mat_info* . ``col`` [ *k* +1]
| |tab| ``if`` ( *mat_info* . ``col`` [ *k* ] == *mat_info* . ``col`` [ *k* +1] )
| |tab| |tab| *mat_info* . ``row`` [ *k* ] < *mat_info* . ``row`` [ *k* +1]

we say that *mat_info* is in column major order.

Row Major Order
===============
If for *k* = 0 , ... , *K* ``-1`` ,

| |tab| *mat_info* . ``row`` [ *k* ] <= *mat_info* . ``row`` [ *k* +1]
| |tab| ``if`` ( *mat_info* . ``row`` [ *k* ] == *mat_info* . ``row`` [ *k* +1] )
| |tab| |tab| *mat_info* . ``col`` [ *k* ] < *mat_info* . ``col`` [ *k* +1]

we say that *mat_info* is in row major order.

Lower Triangular
================
If for *k* = 0 , ... , *K* ``-1`` ,

   *mat_info* . ``row`` [ *k* ] >= *mat_info* . ``col`` [ *k* ]

we say that *mat_info* is lower triangular.

{xrst_end sparse_mat_info}
*/
# include <cppad/utility/vector.hpp>

namespace CppAD { namespace mixed {
   struct sparse_mat_info {
      CppAD::vector<size_t>  row;
      CppAD::vector<size_t>  col;
      CppAD::vector<double>  val;
      //
      void resize(size_t size)
      {  row.resize(size);
         col.resize(size);
         val.resize(size);
      }
   };
} }

# endif
