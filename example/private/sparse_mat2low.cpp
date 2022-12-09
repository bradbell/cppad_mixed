// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin sparse_mat2low.cpp dev}

sparse_mat2low: Example and Test
################################

Private
*******
This example is not part of the
:ref:`cppad_mixed public API<base_class-name>` .

{xrst_literal
   // BEGIN C++
   // END C++
}

{xrst_end sparse_mat2low.cpp}
*/
// BEGIN C++
# include <cstddef>
# include <cppad/mixed/sparse_mat2low.hpp>

bool sparse_mat2low_xam(void)
{  bool ok = true;
   typedef Eigen::SparseMatrix<size_t> sparse_matrix;
   typedef Eigen::Index                Int;
   //
   size_t nr = 3;
   sparse_matrix matrix;
   matrix.resize(Int(nr), Int(nr));
   Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> check(nr, nr);
   size_t count = 0;
   for(size_t i = 0; i < nr; i++)
   {  for(size_t j = 0; j < nr; j++)
      {  matrix.insert(Int(i), Int(j)) = count;
         check(i, j)                   = count;
         ++count;
      }
   }
   sparse_matrix lower = CppAD::mixed::sparse_mat2low(matrix);
   //
   // check the result values
   count = 0;
   for(size_t k = 0; k < nr; k++)
   {  for(sparse_matrix::InnerIterator itr(lower, Int(k)); itr; ++itr)
      {  Int i = itr.row();
         Int j = itr.col();
         ok   &= itr.value() == check(i, j);
         // make sure only lower triangle is inclued
         ok   &= j <= i;
         ++count;
      }
   }
   // make sure all of lower traingle is included
   ok &= count == (nr * (nr + 1) ) / 2;
   return ok;
}
// END C++
