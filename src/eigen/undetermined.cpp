// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin undetermined dev}
{xrst_spell
   nc
   nr
   tol
}

Express An Undetermined Linear System As Dependent and Independent Variables
############################################################################

Syntax
******
*rank* = ``CppAD::mixed::undetermined`` ( *A* , *b* , *tol* , *D* , *I* , *C* , *e* )

Prototype
*********
{xrst_literal
   // BEGIN PROTOTYPE
   // END PROTOTYPE
}

Private
*******
This function is an implementation detail and not part of the
CppAD Mixed user API.

Purpose
*******
We are give a matrix :math:`A \in \B{R}^{m \times n}`
and a vector :math:`b \in \B{R}^m`,
where :math:`m < n` and :math:`A` has rank :math:`m`.
Furthermore, we are interested in the linear constraint equation

.. math::

   A x = b

A matrix :math:`C \in \B{R}^{m \times (n - m)}`
and a vector :math:`e \in \B{R}^m`,
such that the constraint is equivalent to

.. math::

   x_D = C x_I + e

where :math:`D` is a subset, of size :math:`m`,
of the column indices and :math:`I` is the complementary subset of the
column indices.
Note that one can use this equation to
solve for the dependent  variables :math:`x_D`
as a function of independent variables :math:`x_I`.

A
*
It is assumed that the rank of
:math:`A` is equal to *A* . ``rows`` () .

nr
==
We use the notation *nr* = *A* . ``rows`` () ; i.e.,
the number of rows in *A* .
It is ok if *nr*  == 0 no constraints in which case
*D* is empty.

nc
==
We use the notation *nc* = *A* . ``cols`` () ; i.e.,
the number of columns in *A* .

b
*
It is assumed that *b* . ``rows`` () = *nr*
(note that *b* . ``cols`` () == 1 ).

tol
***
This is the tolerance, relative to one,
used for detecting a rank deficient matrix.
To be specific, *tol* times the maximum value in a row
is considered a zero element for that row.

D
*
It is assumed that *D* . ``rows`` () = *nr*
(note that *D* . ``cols`` () == 1 ).
The input value of its elements does not matter.
If *rank* == *nr* ,
upon return the vector :math:`x_D` is

   ( *x* [ *D* [0]] , *x* [ *D* [1]] , ... , *x* [ *D* [ *nr* ``-1`` ]] )^ ``T``

I
*
It is assumed that *I* . ``rows`` () = *nc* ``-`` *nr*
(note that *I* . ``cols`` () == 1 ).
The input value of its elements does not matter.
If *rank* == *nr* ,
upon return the vector :math:`x_I` is

   ( *x* [ *I* [0]] , *x* [ *I* [1]] , ... , *x* [ *I* [ *nr* ``-`` *nc* ``-1`` ]] )^ ``T``

Furthermore the union of the sets corresponding
to :math:`D` and :math:`I` is { 0 , ... , *nc* ``-1`` } .
It follows that the sets do not intersect and none of the elements are
repeated in the vectors *D* or *I* .

C
*
It is assumed that *C* . ``rows`` () == *nr* and
*C* . ``cols`` () == *nc* ``-`` *nr* .
The input value of its elements does not matter.
If *rank* == *nr* ,
upon return it is the matrix :math:`C` in :math:`x_D = C x_I + e`.

e
*
It is assumed that *e* . ``rows`` () == *nr*
(note that *e* . ``cols`` () == 1 ).
The input value of its elements does not matter.
If *rank* == *nr* ,
upon return it is the vector :math:`e` in :math:`x_D = C x_I + e`.

rank
****
The return value has prototype

   ``size_t`` *rank*

and it is the rank of the matrix to tolerance
:ref:`undetermined@tol` .
{xrst_toc_hidden
   example/private/undetermined.cpp
}

Example
*******
The file :ref:`undetermined.cpp-name` is an example
and test of ``undetermined`` .

2DO
***
This routine uses dense matrices, perhaps it would be useful
to convert this (and :ref:`sample_fixed-name` ) to all sparse matrices.

{xrst_end undetermined}
------------------------------------------------------------------------------
*/


# include <cmath>
# include <iostream>
# include <Eigen/Core>

namespace {
   using Eigen::Dynamic;
   typedef Eigen::Matrix<double, Dynamic, Dynamic> double_mat;
   typedef Eigen::Matrix<bool,   Dynamic, 1>       bool_vec;
   typedef Eigen::Matrix<size_t, Dynamic, 1>       size_vec;
   typedef std::pair<size_t, size_t>               size_pair;
   //
   size_pair max_abs(
      const double_mat&    E        ,
      const bool_vec&      row_used ,
      const bool_vec&      col_used )
   {
      // number of rows in E
      size_t nr      = size_t( E.rows() );
      // number of columns in E
      size_t nc      = size_t( E.cols() ) - 1;
      //
      size_pair ret(0,0);
      double max_abs_val = -1.0;
      for(size_t i = 0; i < nr; i++) if( ! row_used[i] )
      {  for(size_t j = 0; j < nc; j++) if( ! col_used[j] )
         {  if( std::fabs( E(i, j) ) > max_abs_val )
            {  ret = size_pair(i, j);
               max_abs_val = std::fabs( E(i, j) );
            }
         }
      }
      // check for case where there is no possible pivot left
      assert( max_abs_val >= 0.0 );
      //
      return ret;
   }
   //
   void elementary(size_pair pivot, double_mat& E)
   {  size_t nr = size_t( E.rows() );
      size_t r  = pivot.first;
      size_t c  = pivot.second;
      double v  = E(r, c);
      E.row(r) /= v;
      // fix roundoff on piovot element
      E(r, c) = 1.0;
      for(size_t i = 0; i < nr; i++)
      {  if( i != r )
         {  E.row(i) -= E(i, c) * E.row(r);
            // fix roundoff on pivot column
         }
      }
   }
}

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN PROTOTYPE
size_t undetermined(
   const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& A     ,
   const Eigen::Matrix<double, Eigen::Dynamic, 1>&              b     ,
   double                                                       tol   ,
   Eigen::Matrix<size_t, Eigen::Dynamic, 1>&                    D     ,
   Eigen::Matrix<size_t, Eigen::Dynamic, 1>&                    I     ,
   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&       C     ,
   Eigen::Matrix<double, Eigen::Dynamic, 1>&                    e     )
// END PROTOTYPE
{  size_t nr = A.rows();
   size_t nc = A.cols();
   assert(  nr < nc );
   assert(  size_t( b.rows() ) == nr );
   assert(  size_t( D.rows() ) == nr );
   assert(  size_t( I.rows() ) == nc - nr );
   assert(  size_t( C.rows() ) == nr );
   assert(  size_t( e.rows() ) == nr );
   assert(  size_t( C.cols() ) == nc - nr );
   //
   // special case
   if( nr == 0 )
   {  for(size_t i = 0; i < nc; i++)
         I[i] = i;
   }
   //
   // E = [ A | b ]
   double_mat E(nr, nc + 1 );
   E.block(0, 0,  nr, nc) = A;
   E.block(0, nc, nr, 1)  = b;
   //
   // Normalize E so all rows have same size maximum element
   for(size_t i = 0; i < nr; i++)
   {  double scale = 0.0;
      for(size_t j = 0; j < nc; j++)
         scale = std::max(scale, fabs( E(i, j) ) );
      if( scale > 0.0 )
         E.row(i) /= scale;
   }
   // note that the maximum absolute element in E(:,nc-1) is 1
   //
   // which rows and colums have been used for pivots
   bool_vec row_used(nr), col_used(nc);
   for(size_t i = 0; i < nr; i++)
      row_used[i] = false;
   for(size_t j = 0; j < nc; j++)
      col_used[j] = false;
   //
   // mapping from pivot row to pivot column
   size_vec pivotrow2col(nr);
   //
   //
   for(size_t rank = 0; rank < nr; rank++)
   {  //
      // determine the next pivot element
      size_pair pivot = max_abs(E, row_used, col_used);
      size_t r = pivot.first;
      size_t c = pivot.second;
      if( std::fabs( E(r, c) ) <= tol )
         return rank;
      //
      // preform elementary row operations for this pivot
      elementary(pivot, E);
      //
      // mark this pivot row and colum as used
      row_used[r] = true;
      col_used[c] = true;
      //
      // record the column corresponding to this row
      pivotrow2col[r] = c;
   }
   //
   // D
   D = pivotrow2col;
   //
   // I and C
   size_t k = 0;
   for(size_t j = 0; j < nc; j++)
   {  // skip columns that are used for pivot operations
      if( ! col_used[j] )
      {  I[k] = j;
         C.col(k) = - E.col(j);
         k++;
      }
   }
   //
   // e
   e = E.col(nc);
   //
   return nr;
}

} } // END_CPPAD_MIXED_NAMESPACE
