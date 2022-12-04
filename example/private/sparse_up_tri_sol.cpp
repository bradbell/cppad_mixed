// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin sparse_up_tri_sol.cpp}

sparse_up_tri_sol: Example and Test
###################################

Private
*******
This example is not part of the
:ref:`cppad_mixed public API<base_class-name>` .

Problem
*******
This example solves the equation

| |tab| [ 1 2 4 ]            [ 1 2 4 ]
| |tab| [ 0 3 5 ] * *result*  = [ 0 3 5 ]
| |tab| [ 0 0 6 ]            [ 0 0 6 ]

{xrst_literal
   // BEGIN C++
   // END C++
}

{xrst_end sparse_up_tri_sol.cpp}
*/
// BEGIN C++
# include <cppad/mixed/sparse_up_tri_sol.hpp>
# include <iostream>

bool sparse_up_tri_sol_xam(void)
{  bool ok = true;
   double eps = 100. * std::numeric_limits<double>::epsilon();
   //
   typedef Eigen::SparseMatrix<double, Eigen::RowMajor> sparse_by_row;
   typedef Eigen::SparseMatrix<double, Eigen::ColMajor> sparse_by_col;

   // left
   size_t nr = 3;
   sparse_by_row left;
   left.resize(int(nr), int(nr));
   size_t count = 0;
   for(size_t i = 0; i < nr; i++)
   {  for(size_t j = i; j < nr ; j++)
         left.insert(int(i), int(j)) = double( ++count );
   }
   // right
   size_t nc = 3;
   sparse_by_col right;
   right.resize(int(nr), int(nc));
   right = left;
   //
   sparse_by_col result = CppAD::mixed::sparse_up_tri_sol(left, right);
   //
   // check that the result is a upper triangular verison of identity matrix
   for(size_t j = 0; j < nc; j++)
   {  count = 0;
      for(sparse_by_col::InnerIterator itr(result, int(j)); itr; ++itr)
      {  ++count;
         double check = 0.0;
         if( itr.row() == itr.col() )
            check = 1.0;
         ok &= std::fabs( itr.value() - check ) <= eps;
         ok &= itr.row() <= itr.col();
      }
   }
   return ok;
}
// END C++
