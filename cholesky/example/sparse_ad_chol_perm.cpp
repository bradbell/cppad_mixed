// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include "../sparse_ad_cholesky.hpp"

/*
{xrst_begin sparse_ad_chol_perm.cpp dev}

Sparse AD Cholesky Permutation: Example and Test
################################################

Source
******
{xrst_literal
   // BEGIN C++
   // END C++
}

{xrst_end sparse_ad_chol_perm.cpp}
*/
// BEGIN C++
bool sparse_ad_chol_perm(void)
{  using CppAD::AD;
   //
   bool ok        = true;
   AD<double> eps = 100. * std::numeric_limits<double>::epsilon();
   // --------------------------------------------------------------------
   // create sparse_ad_cholesky object
   int nc = 3;
   Eigen::SparseMatrix< AD<double>, Eigen::ColMajor> ad_Blow(nc, nc);
   ad_Blow.insert(0,0) = 1.0; //     [ 1.0   0.0    0.5 ]
   ad_Blow.insert(2,0) = 0.5; // B = [ 0.0   0.5    0.0 ]
   ad_Blow.insert(1,1) = 0.5; //     [ 0.5   0.0    2.0 ]
   ad_Blow.insert(2,2) = 2.0;
   CppAD::mixed::sparse_ad_cholesky cholesky;
   cholesky.initialize( ad_Blow );
   //
   // Permutation matgrix kcorresponding to this cholesky
   const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>& P =
      cholesky.permutation();
   //
   // Lower triangle of symmetric matrix with same sparsity pattern as B
   Eigen::SparseMatrix< AD<double>, Eigen::ColMajor> Alow(nc, nc);
   Alow.insert(0,0) = 2.0; //     [ 2.0   0.0   0.5  ]
   Alow.insert(2,0) = 0.5; // A = [ 0.0   0.5   0.0  ]
   Alow.insert(1,1) = 0.5; //     [ 0.5   0.0   1.0  ]
   Alow.insert(2,2) = 1.0;
   //
   // compute the Choleksy factorization of A
   Eigen::SparseMatrix< AD<double>, Eigen::ColMajor> L;
   cholesky.eval(Alow, L);
   ok &= L.rows() == nc;
   ok &= L.cols() == nc;
   //
   // compute P^T * L * L^T * P
   Eigen::SparseMatrix< AD<double>, Eigen::ColMajor> tmp1, tmp2;
   tmp1 = P.transpose() * L;
   tmp2 = tmp1 * L.transpose();
   tmp1 = tmp2 * P;
   //
   // check that A = P^T * L * L^T * P
   Eigen::Matrix< AD<double>, 3, 3> A( Alow ), prod( tmp1 );
   for(size_t i = 0; i < 3; i++)
   {  for(size_t j = 0; j < 3; j++)
      {  if( j > i )
            A(i, j) = A(j, i);
         ok &= CppAD::NearEqual( A(i, j), prod(i, j), eps, eps );
      }
   }
   // -----------------------------------------------------------------------
   return ok;
}
// END C++
