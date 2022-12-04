// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin sparse_low_tri_sol}
{xrst_spell
   cols
}

Solve a Sparse Lower Triangular Linear System
#############################################

Syntax
******
*result* = ``CppAD::mixed::sparse_low_tri_sol`` ( *left* , *right* ) .

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

left
****
This must be a square invertible lower triangular matrix; i.e.,
*left* . ``rows`` () == *left* . ``cols`` () and
for each *left* ( *i* , *j* ) in this sparse matrix,
*i* >= *j* .

right
*****
The number of rows in this matrix must equal the number of columns in
*left* ; i.e.,
*right* . ``rows`` () == *left* . ``cols`` () .

result
******
This result has the same number of rows as *left* ,
the same number of columns as *right* ,
and satisfies the equation

   *left* * *result* = *right*

where ``*`` is matrix multiplication.
{xrst_toc_hidden
   example/private/sparse_low_tri_sol.cpp
}
Example
*******
The file :ref:`sparse_low_tri_sol.cpp-name` is an example
and test of ``sparse_low_tri_sol`` .

{xrst_end sparse_low_tri_sol}
*/
# include <cppad/mixed/sparse_low_tri_sol.hpp>
# include <iostream>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

/*
template <class sparse_type>
void print(const std::string& label, const sparse_type& mat)
{  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m = mat;
   std::cout << label << "=\n" << m << "\n";
}
*/


// BEGIN PROTOTYPE
Eigen::SparseMatrix<double, Eigen::ColMajor> sparse_low_tri_sol(
   const Eigen::SparseMatrix<double, Eigen::RowMajor>&  left  ,
   const Eigen::SparseMatrix<double, Eigen::ColMajor>&  right )
// END PROTOTYPE
{  using Eigen::ColMajor;
   using Eigen::RowMajor;
   using Eigen::Dynamic;
   typedef Eigen::Index Int;
   //
   typedef Eigen::SparseMatrix<double, ColMajor>::InnerIterator col_itr;
   typedef Eigen::SparseMatrix<double, RowMajor>::InnerIterator row_itr;
   //
   // number of rows and columns in the result
   Int nr = left.rows();
   Int nc = right.cols();
   assert( nr == left.cols() );
   assert( nr == right.rows() );
   //
   // initialize result matrix as empty
   Eigen::SparseMatrix<double, ColMajor> result(nr, nc);
   //
   // One column of the result as a singly linked list with no
   // memory allocation when adding and removing elements.
   // initialize it as zero
   Int res_first = nr; // null pointer value
   Eigen::Matrix<double, Dynamic, 1> res_value(nr);
   Eigen::Matrix<Int, Dynamic, 1>    res_next(nr);
   res_value.setZero();
   //
   // One column of the right hand side as a singly linked list
   // with no memory allocation when adding and removing elements.
   Int rhs_first;
   Eigen::Matrix<double, Dynamic, 1> rhs_value(nr);
   Eigen::Matrix<Int, Dynamic, 1>    rhs_next(nr);
   //
   // for each column in the right matrix
   for(Int j = 0; j < right.cols(); ++j)
   {  // iterator for non-zero entries in j-th column of right hand side
      col_itr right_itr(right, j);
      //
      // -----------------------------------------------------------------
      // convert this column of the right hand side to a singly linked list
      rhs_first = nr;   // null pointer
      if( right_itr )
      {  assert( right_itr.col() == j );
         //
         Int i        = right_itr.row();
         rhs_value[i] = right_itr.value();
         rhs_first    = i;
         while( ++right_itr )
         {  assert( right_itr.col() == j );
            //
            rhs_next[i]  = right_itr.row();
            i            = right_itr.row();
            rhs_value[i] = right_itr.value();
         }
         rhs_next[i] = nr; // null pointer
      }
      // -----------------------------------------------------------------
      //
      // initialize row index for right hand side that is greater than
      // or equal the current row being solved
      Int rhs_ge = rhs_first;
      //
      // initialize pointer to previous non-zero entry in result
      Int res_previous = nr; // null pointer
      //
      // for each row in the left matrix
      for(Int i = rhs_first; i < nr; ++i)
      {  // advance to next right hand side index greater than or equal i
         while( rhs_ge < i )
            rhs_ge = rhs_next[ rhs_ge ];
         //
         // initialize summation for (i, j) entry of the result
         if( rhs_ge == i )
         {  // --------------------------------------------------------
            // this means we have a non-zero result at index i
            if( res_previous == nr )
               res_first = i;
            else
               res_next[res_previous] = i;
            //
            res_value[i] = rhs_value[i];
            res_next[i]  = nr;
            res_previous = i;
            // --------------------------------------------------------
         }
         //
         // row rhs_first must have a non-zero result
         assert( res_previous != nr );
         //
         // for each entry in the i-th row of left matrix
         double left_ii = 0.0;
         for(row_itr left_itr(left, i); left_itr; ++left_itr)
         {  // (i, k) index in left matrix
            Int k = left_itr.col();
            // check that left is lower triangular
            assert( i >= k );
            //
            if( i > k )
            {  if (res_value[k] != 0.0 )
               {  // ----------------------------------------------------
                  // this means that we have a non-zero result at index i
                  res_value[i] -= left_itr.value() * res_value[k];
                  if( res_previous != i )
                  {  res_next[res_previous] = i;
                     res_next[i]            = nr;
                     res_previous           = i;
                  }
                  // ----------------------------------------------------
               }
            }
            else
               left_ii = left_itr.value();
         }
         assert( left_ii != 0.0 );
         //
         // check if we have a non-zero entry to divide
         if( res_previous == i )
         {  res_value[i] /= left_ii;
            result.insert(i, j) = res_value[i];
         }
      }
      // restrore the res_value vector to all zeros
      Int i = res_first;
      while( i < nr )
      {  res_value[i] = 0.0;
         i            = res_next[i];
      }
      // restor the res list to empty
      res_first = nr;
   }
   return result;
}
} } // END_CPPAD_MIXED_NAMESPACE
