// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin sparse_scale_diag.cpp}
{xrst_spell
   diag
}

sparse_scale_diag: Example and Test
###################################

Private
*******
This example is not part of the
:ref:`cppad_mixed public API<base_class-name>` .

Description
***********
The call to :ref:`sparse_scale_diag-name` below computes the matrix

| |tab| [ 3* 1 0   0   ]
| |tab| [ 2   3* 3 0   ]
| |tab| [ 4   5   3* 6 ]

{xrst_literal
   // BEGIN C++
   // END C++
}

{xrst_end sparse_scale_diag.cpp}
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
