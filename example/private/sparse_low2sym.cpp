// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin sparse_low2sym.cpp}
{xrst_spell
   sym
}

sparse_low2sym: Example and Test
################################

Private
*******
This example is not part of the
:ref:`cppad_mixed public API<base_class-name>` .

{xrst_literal
   // BEGIN C++
   // END C++
}

{xrst_end sparse_low2sym.cpp}
*/
// BEGIN C++
# include <cmath>
# include <cstddef>
# include <limits>
# include <cppad/mixed/sparse_low2sym.hpp>

bool sparse_low2sym_xam(void)
{  bool ok = true;
   typedef Eigen::SparseMatrix<size_t> sparse_matrix;
   typedef Eigen::Index                Int;
   //
   size_t nr = 3;
   sparse_matrix lower;
   lower.resize(Int(nr), Int(nr));
   Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> check(nr, nr);
   size_t count = 0;
   for(size_t i = 0; i < nr; i++)
   {  for(size_t j = 0; j <= i; j++)
      {  lower.insert(Int(i), Int(j)) = count;
         check(i, j)                  = count;
         check(j, i)                  = count;
         ++count;
      }
   }
   sparse_matrix symmetric = CppAD::mixed::sparse_low2sym(lower);
   //
   // check the result
   count = 0;
   for(size_t k = 0; k < nr; k++)
   {  for(sparse_matrix::InnerIterator itr(symmetric, Int(k)); itr; ++itr)
      {  Int i = itr.row();
         Int j = itr.col();
         ok   &= itr.value() == check(i, j);
         ++count;
      }
   }
   ok &= count == nr * nr;
   return ok;
}
// END C++
