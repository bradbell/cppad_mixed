// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
A      = [ 5, 4, 2; 4, 5, 1; 2, 1, 5]
det(A) = 36
inv(A) = [ 24, -18, -6; -18, 21, 3; -6, 3, 9] / det(A)
*/
# include <Eigen/Sparse>
# include <vector>
# include <iostream>
# include <cppad/cppad.hpp>

namespace { // BEGIN empty namespace

bool test_one(void)
{  bool ok = true;
   double eps = 100. * std::numeric_limits<double>::epsilon();

   using Eigen::Dynamic;
   typedef CppAD::AD<double>                     real;
   typedef Eigen::SparseMatrix<real>             real_sparse_matrix;
   typedef Eigen::Matrix<real, Dynamic, Dynamic> real_dense_matrix;

   // size of the matrix
   size_t n = 3;

   double A_inv[] = {
       24.0, -18.0, -6.0,
      -18.0,  21.0,  3.0,
       -6.0,   3.0,  9.0
   };
   for(size_t i = 0; i < sizeof(A_inv)/sizeof(A_inv[0]); i++)
      A_inv[i] /= 36.;

   // Start tape of function that maps non-zeros in lower triangle of a
   // positive definate A to the log of its determinant
   CppAD::vector<real> A_vec(n * (n + 1) / 2 );
   A_vec[0] = 5.0; // A(0, 0)
   A_vec[1] = 4.0; // A(1, 0)
   A_vec[2] = 5.0; // A(1, 1)
   A_vec[3] = 2.0; // A(2, 0)
   A_vec[4] = 1.0; // A(2, 1)
   A_vec[5] = 5.0; // A(2, 2)
   CppAD::Independent( A_vec );


   // convert to an Eigen sparse matrix
   real_sparse_matrix A;
   A.resize(int(n), int(n));

   // create Matrix
   A.insert(0, 0) = A_vec[0];
   A.insert(1, 0) = A_vec[1];
   A.insert(1, 1) = A_vec[2];
   A.insert(2, 0) = A_vec[3];
   A.insert(2, 1) = A_vec[4];
   A.insert(2, 2) = A_vec[5];

   // Compute the log of determinant using lower triangle of A
   Eigen::SimplicialLDLT<real_sparse_matrix, Eigen::Lower> chol;
   chol.analyzePattern(A);
   chol.factorize(A);
   real_dense_matrix diag = chol.vectorD();
   ok &= size_t( diag.size() ) == n;
   CppAD::vector<real> log_det(1);
   log_det[0] = 0.0;
   for(size_t i = 0; i < n; i++)
      log_det[0] += log( diag(i) );

   // Stop tape and create function
   CppAD::ADFun<double>(A_vec, log_det);

   // check function value
   ok &= fabs( exp(log_det[0]) / 36. - 1.0 ) < eps;

   // initialize dense column vector as zero
   real_dense_matrix b(n, 1);
   for(size_t i = 0; i < n; i++)
      b(i) = 0.0;

   // check computation of the inverse of A
   for(size_t j = 0; j < n; j++)
   {  b(j) = 1.0;
      real_dense_matrix x = chol.solve(b);
      b(j) = 0.0;
      for(size_t i = 0; i < n; i++)
         ok &= fabs( x(i) / A_inv[i*n + j] - 1.0 ) <= eps;
   }
   return ok;
}

bool test_two(void)
{  bool ok = true;
   double eps = 100. * std::numeric_limits<double>::epsilon();

   using Eigen::Dynamic;
   typedef double                                   real;
   typedef Eigen::SparseMatrix<real>                real_sparse_matrix;
   typedef Eigen::SparseMatrix<real>::InnerIterator inner_itr;

   // size of the matrix
   size_t n = 3;

   double A_inv[] = {
       24.0, -18.0, -6.0,
      -18.0,  21.0,  3.0,
       -6.0,   3.0,  9.0
   };
   for(size_t i = 0; i < sizeof(A_inv)/sizeof(A_inv[0]); i++)
      A_inv[i] /= 36.;

   CppAD::vector<real> A_vec(n * (n + 1) / 2 );
   A_vec[0] = 5.0; // A(0, 0)
   A_vec[1] = 4.0; // A(1, 0)
   A_vec[2] = 5.0; // A(1, 1)
   A_vec[3] = 2.0; // A(2, 0)
   A_vec[4] = 1.0; // A(2, 1)
   A_vec[5] = 5.0; // A(2, 2)

   // Lower triangle of positive definate A
   real_sparse_matrix A;
   A.resize(int(n), int(n));
   A.insert(0, 0) = 5.0;
   A.insert(1, 0) = 4.0;
   A.insert(1, 1) = 5.0;
   A.insert(2, 0) = 2.0;
   A.insert(2, 1) = 1.0;
   A.insert(2, 2) = 5.0;

   // Compute the log of determinant using lower triangle of A
   Eigen::SimplicialLDLT<real_sparse_matrix, Eigen::Lower> chol;
   chol.analyzePattern(A);
   chol.factorize(A);

   // check computation of the inverse of A
   size_t count = 0;
   for(size_t j = 0; j < n; j++)
   {  real_sparse_matrix b;
      b.resize( int(n), 1);
      b.insert( int(j), 0) = 1.0;
      real_sparse_matrix x = chol.solve(b);
      for(int k = 0; k < x.outerSize(); ++k)
      {  ok &= k == 0;
         for(inner_itr itr(x, k); itr; ++itr)
         {  size_t i = itr.row();
            real value = itr.value();
            ok &= std::fabs( value - A_inv[ i * n + j] ) < 10. * eps;
            count++;
         }
      }
   }
   ok &= count == 9;

   return ok;
}

} // END namespace

bool eigen_xam(void)
{  bool ok  = true;
   ok &= test_one();
   ok &= test_two();
   return ok;
}
