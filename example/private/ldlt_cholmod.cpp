// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-25 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin ldlt_cholmod.cpp dev}
{xrst_spell
  ccc
  cov
  logdet
  nrow
  sim
  rcond
}

Example Using ldlt_cholmod Class
################################

Problem Description
*******************
We define the lower triangular matrix

.. math::

   L =
   \left( \begin{array}{ccc}
      1 & 0 & 0 \\
      1 & 1 & 0 \\
      1 & 1 & 1
   \end{array} \right)
   \W{,}
   D =
   \left( \begin{array}{ccc}
      1 & 0 & 0 \\
      0 & 4 & 0 \\
      0 & 0 & 9
   \end{array} \right)
   \W{,}
   L^\R{T} =
   \left( \begin{array}{ccc}
      1 & 1 & 1 \\
      0 & 1 & 1 \\
      0 & 0 & 1
   \end{array} \right)

and the positive definite matrix

.. math::

   H = L D L^\R{T} =
   \left( \begin{array}{ccc}
      1 & 1 & 1 \\
      1 & 5 & 5 \\
      1 & 5 & 14
   \end{array} \right)

The inverse of :math:`H` is given by

.. math::

   H^{-1} = L^\R{-T} D^{-1} L^{-1} =
   \frac{1}{36}
   \left( \begin{array}{ccc}
      45  & -9  & 0  \\
      -9  & 13  & -4 \\
      0   & -4  & 4
   \end{array} \right)

which can be checked by multiplying by :math:`H H^{-1}`.

constructor
***********
See the following code below:
::

   CppAD::mixed::ldlt_cholmod ldlt_obj(nrow);

init
****
See the following under
:ref:`ldlt_cholmod.cpp@Source Code` below:
::

   ldlt_obj.init(H_rcv.pat());

pattern
*******
See the following under
:ref:`ldlt_eigen.cpp@Source Code` below:
::

   H_rc = ldlt_obj.pattern();

update
******
See the following under Source Code below:
::

   ldlt_obj.update(H_rcv);

rcond
*****
See the following under Source Code below:
::

   rcond_D = ldlt_obj.rcond( );

logdet
******
See the following under Source Code below:
::

   logdet_H = ldlt_obj.logdet(negative);

solve_H
*******
See the following under Source Code below:
::

   ldlt_obj.solve_H(row, val_in, val_out);

sim_cov
*******
See the following under Source Code below:
::

   ok &= ldlt_obj.sim_cov(w, v)

Source Code
***********
{xrst_literal
   // BEGIN C++
   // END C++
}

{xrst_end ldlt_cholmod.cpp}
*/
// BEGIN C++
# include <cppad/mixed/ldlt_cholmod.hpp>
# include <limits>
# include <cmath>
# include <cassert>

bool ldlt_cholmod_xam(void)
{  bool ok    = true;
   double eps = 100. * std::numeric_limits<double>::epsilon();

   double H_inv[] = {
      45.0,  -9.0,  0.0,
      -9.0,  13.0, -4.0,
      0.0,   -4.0,  4.0
   };
   for(size_t i = 0; i < sizeof(H_inv)/sizeof(H_inv[0]); i++)
      H_inv[i] /= 36.;

   // create cholmod object
   size_t nrow = 3;    // number of rows in H
   size_t ncol = nrow; // number of columns in H
   CppAD::mixed::ldlt_cholmod ldlt_obj(nrow);
   assert( nrow * ncol == sizeof(H_inv) / sizeof(H_inv[0]) );

   // sparsity pattern for the lower triangular of H (dense)
   size_t nnz = 6;
   CppAD::mixed::sparse_rc H_rc(nrow, ncol, nnz);
   {  size_t k = 0;
      H_rc.set(k++, 0, 0);
      H_rc.set(k++, 1, 0);
      H_rc.set(k++, 2, 0);
      H_rc.set(k++, 1, 1);
      H_rc.set(k++, 2, 1);
      H_rc.set(k++, 2, 2);
   }

   // values in lower triangle of H
   CppAD::mixed::d_sparse_rcv H_rcv( H_rc );
   {  size_t k = 0;
      H_rcv.set(k++, 1.0);  // H_0,0  = 1.0
      H_rcv.set(k++, 1.0);  // H_1,0  = 1.0
      H_rcv.set(k++, 1.0);  // H_2,0  = 1.0
      H_rcv.set(k++, 5.0);  // H_1,1  = 5.0
      H_rcv.set(k++, 5.0);  // H_2,1  = 5.0
      H_rcv.set(k++, 14.0); // H_2,2  = 14.0
   }
   //
   // initialize the matrix using only the sparsity pattern
   ldlt_obj.init( H_rcv.pat() );

   // check the pattern function
   const CppAD::mixed::sparse_rc& pattern( ldlt_obj.pattern() );
   ok &= pattern.nnz() == nnz;
   for(size_t k = 0; k < nnz; ++k)
   {  ok &= pattern.row()[k] == H_rc.row()[k];
      ok &= pattern.col()[k] == H_rc.col()[k];
   }

   // factor the matrix using the values
   ldlt_obj.update( H_rcv );

   // compute log of determinant of H
   size_t negative;
   double logdet_H = ldlt_obj.logdet(negative);
   ok &= negative == 0;
   ok &= std::fabs( logdet_H / std::log(36.0) - 1.0 ) <= eps;

   // ok
   // This test will fail if the Permutation is not the identity;
   // see ldlt_rcond in the documentation.
   double rcond_D = ldlt_obj.rcond( );
   ok &= std::fabs( rcond_D * 9.0 - 1.0 ) <= eps;

   // test solve_H
   {  CppAD::vector<size_t> row(3);
      CppAD::vector<double> val_in(3), val_out(3);
      for(size_t j = 0; j < ncol; j++)
      {  // solve for the j-th column of the inverse matrix
         for(size_t k = 0; k < 3; k++)
         {  row[k] = k;
            if( row[k] == j )
               val_in[k] = 1.0;
            else
               val_in[k] = 0.0;
         }
         ldlt_obj.solve_H(row, val_in, val_out);
         //
         for(size_t k = 0; k < row.size(); k++)
         {  size_t i       = row[k];
            double check_i = H_inv[ i * nrow + j ];
            ok &= std::fabs( val_out[k] - check_i ) <= eps;
         }
      }
   }

   // test inv
   {  size_t K = 6;
      CppAD::vector<size_t> row_in(K), col_in(K);
      CppAD::vector<double> val_out(K);
      {  size_t k = 0;
         // use row major ordering to test sorting
         for(size_t i = 0; i < nrow; i++)
         {  // only request lower triangle
            for(size_t j = 0; j <= i; j++)
            {  row_in[k] = i;
               col_in[k] = j;
               k++;
            }
         }
      }
      ldlt_obj.inv(row_in, col_in, val_out);
      for(size_t k = 0; k < K; k++)
      {  double check = H_inv[ row_in[k] * nrow + col_in[k] ];
         ok          &= std::fabs( val_out[k] - check ) <= eps;
      }
   }
   //
   // test sim_cov
   CppAD::vector<double> w(3), v(3), c(3);
   for(size_t i = 0; i < 3; i++)
      w[i] = double( 2 * i + 1);
   ok &= ldlt_obj.sim_cov(w, v);
   // solve w = L^{T} c
   c[2] = w[2] / 3.0;
   c[1] = ( w[1] - 2 * c[2] ) / 2.0;
   c[0] = ( w[0] - 1.0 * c[1] - 1.0 * c[2] ) / 1.0;
   for(size_t i = 0; i < 3; i++)
      ok  &= std::fabs( v[i] - c[i] ) <= eps;

   return ok;
}
// END C++
