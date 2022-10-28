// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/typedef.hpp>
# include "../sparse_ad_cholesky.hpp"

/*
$begin sparse_ad_chol_var.cpp$$
$spell
   Cholesky
$$

$section Sparse AD Cholesky Variable Calculation: Example and Test$$

$head Problem$$
We are given the function $latex A : \B{R}^3 \rightarrow \B{R}^{3 \times 3}$$
defined by
$latex \[
   A(x) = \left( \begin{array}{ccc}
      x_0 & 0    & x_2  \\
      0   & x_1  & 0   \\
      x_2 & 0    & x_3
   \end{array} \right)
\] $$

$head Permutation$$
The fill reducing permutation
$cref/P/sparse_ad_cholesky/Notation/P/$$
transposes indices zero and one; i.e.
$latex \[
   P= \left( \begin{array}{ccc}
      0   & 1    & 0    \\
      1   & 0    & 0   \\
      0   & 0    & 1
   \end{array} \right)
   \; , \;
   P A(x) P^\R{T} = \left( \begin{array}{ccc}
      x_1 & 0    & 0    \\
      0   & x_0  & x_2 \\
      0   & x_2  & x_3
   \end{array} \right)
\] $$

$head Cholesky Factor$$
The Cholesky factor
$cref/L/sparse_ad_cholesky/Notation/L/$$ is
$latex \[
   L(x) = \left( \begin{array}{ccc}
      \sqrt{x_1} & 0                   & 0    \\
      0          & \sqrt{x_0}          & 0   \\
      0          & x_2 / \sqrt{x_0}    & \sqrt{ x_3 - x_2^2 / x_0 }
   \end{array} \right)
\] $$
This can be verified by checking
$latex P * A(x) * P^\R{T} = L * L^\R{T}$$.

$head Source$$
$srcthisfile%4%// BEGIN C++%// END C++%1%$$
$end
*/
// BEGIN C++
bool sparse_ad_chol_var(void)
{  using CppAD::AD;
   using CppAD::mixed::d_vector;
   using CppAD::mixed::a1_vector;
   typedef Eigen::SparseMatrix< AD<double>, Eigen::ColMajor> sparse_ad_matrix;
   //
   bool ok     = true;
   // --------------------------------------------------------------------
   size_t nx = 4;
   size_t ny = 4;
   d_vector x(nx), y(ny);
   a1_vector ax(nx);
   ax[0] = x[0] = 2.0;
   ax[1] = x[1] = 1.0;
   ax[2] = x[2] = 0.5;
   ax[3] = x[3] = 3.0;
   // --------------------------------------------------------------------
   // create sparse_ad_cholesky object
   size_t nc = 3;
   Eigen::SparseMatrix< AD<double>, Eigen::ColMajor> ad_Blow;
   ad_Blow.resize(int(nc), int(nc));
   ad_Blow.insert(0,0) = 2.0;
   ad_Blow.insert(2,0) = 1.0;
   ad_Blow.insert(1,1) = 0.2;
   ad_Blow.insert(2,2) = 3.0;
   CppAD::mixed::sparse_ad_cholesky cholesky;
   cholesky.initialize( ad_Blow );
   // ----------------------------------------------------------------------
   for(size_t k = 0; k < nx; k++)
   {  // Create function object corresponding to f(x)
      CppAD::Independent( ax );
      //
      // vector with only k-th component a variable
      a1_vector ac(nx);
      for(size_t j = 0; j < nx; j++)
      {  if( j == k )
            ac[j] = ax[j];
         else
            ac[j] = x[j];
      }
      //
      // Lower triangle of symmetric matrix with same sparsity pattern as B
      sparse_ad_matrix Alow;
      Alow.resize(int(nc), int(nc));
      Alow.insert(0,0) = ac[0];
      Alow.insert(1,1) = ac[1];
      Alow.insert(2,0) = ac[2];
      Alow.insert(2,2) = ac[3];
      //
      // compute the Choleksy factorization of A
      Eigen::SparseMatrix< AD<double>, Eigen::ColMajor> L;
      cholesky.eval(Alow, L);
      ok &= size_t(L.rows()) == nc;
      ok &= size_t(L.cols()) == nc;
      //
      // check which elements of L are variables
      size_t count = 0;
      for(size_t j = 0; j < nc; j++)
      {  for(sparse_ad_matrix::InnerIterator itr(L, int(j)); itr; ++itr)
         {  ++count;
            size_t i = size_t( itr.row() );
            assert( size_t(itr.col()) == j );
            if( i == 0 && j == 0 )
            {  if( k == 1 )
                  ok &= Variable( itr.value() );
               else
                  ok &= Parameter( itr.value() );
            }
            else if( i == 1 && j == 1 )
            {  if( k == 0 )
                  ok &= Variable( itr.value() );
               else
                  ok &= Parameter( itr.value() );
            }
            else if( i == 2 && j == 1 )
            {  if( k == 0 || k == 2)
                  ok &= Variable( itr.value() );
               else
                  ok &= Parameter( itr.value() );
            }
            else if( i == 2 && j == 2 )
            {  if( k == 0 || k == 2 || k == 3 )
                  ok &= Variable( itr.value() );
               else
                  ok &= Parameter( itr.value() );
            }
            else
               ok = false;
         }
      }
      ok &= count == 4;
      //
      // abort the recording started by k-th call to Independent
      AD<double>::abort_recording();
   }
   // -----------------------------------------------------------------------
   return ok;
}
// END C++
