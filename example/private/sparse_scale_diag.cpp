// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
$begin sparse_scale_diag.cpp$$
$spell
   tri
   cppad
$$

$section sparse_scale_diag: Example and Test$$

$head Private$$
This example is not part of the
$cref/cppad_mixed public API/base_class/$$.

$head Description$$
The call to $cref sparse_scale_diag$$ below computes the matrix
$codei%
   [ 3*1 0   0   ]
   [ 2   3*3 0   ]
   [ 4   5   3*6 ]
%$$

$code
$srcthisfile%0%// BEGIN C++%// END C++%1%$$
$$
$end
*/
// BEGIN C++
# include <cmath>
# include <cstddef>
# include <limits>
# include <cppad/mixed/sparse_scale_diag.hpp>

bool sparse_scale_diag_xam(void)
{  bool ok = true;
   double eps = 10. * std::numeric_limits<double>::epsilon();
   typedef Eigen::SparseMatrix<double> sparse_matrix;
   typedef Eigen::Index                Int;
   //
   Int nr = 3;
   sparse_matrix matrix(nr, nr);
   Int count = 0;
   for(Int i = 0; i < nr; i++)
   {  for(Int j = 0; j <= i; j++)
         matrix.insert(i, j) = double( ++count );
   }
   Int scale = 3;
   CppAD::mixed::sparse_scale_diag(scale, matrix);
   //
   // Check that the diagonal has been scaled by there and that the
   // matrx is lower triangular
   for(Int k = 0; k < nr; k++)
   {  for(sparse_matrix::InnerIterator itr(matrix, k); itr; ++itr)
      {  Int i = itr.row();
         Int j = itr.col();
         ok   &= j <= i;
         if( i == j )
         {  Int original_value = (i + 1) * (i + 2) / 2;
            Int new_value      = scale * original_value;
            ok &= std::fabs( itr.value() - double(new_value) ) < eps;
         }
      }
   }
   return ok;
}
// END C++
