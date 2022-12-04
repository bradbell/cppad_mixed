// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin sparse_eigen2info.cpp}

sparse_eigen2info: Example and Test
###################################

Private
*******
This example is not part of the
:ref:`cppad_mixed public API<base_class-name>` .

{xrst_literal
   // BEGIN C++
   // END C++
}

{xrst_end sparse_eigen2info.cpp}
*/
// BEGIN C++
# include <cppad/mixed/sparse_eigen2info.hpp>

bool sparse_eigen2info_xam(void)
{  bool ok = true;
   typedef Eigen::SparseMatrix<double, Eigen::ColMajor, int> sparse_matrix;
   //
   int nr = 3;
   int nc = nr;
   sparse_matrix matrix(nr, nr);
   size_t count = 0;
   for(int j = 0; j < nc; j++)
   {  for(int i = j; i < nr; i++)
      {  matrix.insert(i, j) = double(count);
         ++count;
      }
   }
   CppAD::mixed::sparse_mat_info info; // empty matrix
   CppAD::mixed::sparse_eigen2info(matrix, info);
   //
   // check the result values
   count = 0;
   for(size_t j = 0; j < size_t(nc); j++)
   {  for(size_t i = j; i < size_t(nr); i++)
      {  ok &= info.row[count] == i;
         ok &= info.col[count] == j;
         ok &= info.val[count] == double(count);
         ++count;
      }
   }
   ok &= count == info.row.size();
   ok &= count == info.col.size();
   ok &= count == info.val.size();
   //
   return ok;
}
// END C++
