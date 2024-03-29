// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin cholmod_factor_xam dev}

Example Using Cholmod LDLT Factors
##################################

Problem Description
*******************
Solve for and check
:math:`L`, :math:`D`, and :math:`P`
in the factorization

.. math::

   P * A * P^{\R{T}} = L * D * L^{\R{T}}

where :math:`A` is a symmetric matrix,
:math:`L` is lower triangular with ones on the diagonal,
:math:`D` is diagonal,
and :math:`P` is a permutation matrix.
{xrst_literal
   // BEGIN C++
   // END C++
}

{xrst_end cholmod_factor_xam}
*/
// BEGIN C++
# include <cppad/mixed/include_cholmod.hpp>
# include <limits>
# include <cmath>
# include <cassert>
# include <vector>

# define CHOLMOD_TRUE                  1
# define CHOLMOD_FALSE                 0
# define CHOLMOD_STYPE_LOWER_TRIANGLE -1


namespace { // BEGIN_EMPTY_NAMESPACE
   void add_T_entry(cholmod_triplet *T, int r, int c, double x)
   {  size_t k           = T->nnz;
      ((int*)T->i)[k]    = r;
      ((int*)T->j)[k]    = c;
      ((double*)T->x)[k] = x;
      T->nnz++;
   }
   void matrix_prod(
      size_t n,
      bool A_transpose,
      const double* A,
      bool B_transpose,
      const double* B,
      double *C
   )
   {  for(size_t i = 0; i < n; ++i)
      {  for(size_t j = 0; j < n; ++j)
         {  C[i * n + j] = 0;
            for(size_t k = 0; k < n; ++k)
            {  double aik = A[i * n + k];
               if( A_transpose )
                  aik = A[k * n + i];
               //
               double bkj = B[k * n + j];
               if( B_transpose )
                  bkj = B[j * n + k];
               //
               C[i * n + j] += aik * bkj;
            }
         }
      }
      return;
   }
}

bool cholmod_factor_xam(void)
{  bool ok      = true;
   //
   // The matrix A
   size_t n = 4;
   const double A[] = {
      1.0, 0.0, 0.0,  1.0,
      0.0, 1.0, 0.0,  0.0,
      0.0, 0.0, 1.0,  0.0,
      1.0, 0.0, 0.0, 11.0
   };
   //
   // initialize common
   cholmod_common com;
   cholmod_start(&com);
   // always do simplicial factorization
   com.supernodal = CHOLMOD_SIMPLICIAL;
   // do LDLT factorization and leave in LDLT form
   com.final_ll = CHOLMOD_FALSE;

   // allocate triplet
   size_t nrow = n;
   size_t ncol = n;
   size_t nzmax = nrow * ncol;
   int    T_stype = CHOLMOD_STYPE_LOWER_TRIANGLE;
   int    T_xtype = CHOLMOD_REAL;
   cholmod_triplet* triplet =
      cholmod_allocate_triplet(nrow, ncol, nzmax, T_stype, T_xtype, &com);
   ok &= triplet->nnz ==  0;

   // triplet entries corresponding to lower triangle of H
   for(size_t i = 0; i < nrow; i++)
   {  for(size_t j = 0; j <= i; j++)
      if( A[i * n + j] != 0.0 )
         add_T_entry(triplet, int(i), int(j), A[i * n + j] );
   }
   //
   // convert triplet to sparse representation of H
   cholmod_sparse* sparse = cholmod_triplet_to_sparse(triplet, 0, &com);
   //
   //  the matrix
   cholmod_factor* factor = cholmod_analyze(sparse, &com);
   int flag               = cholmod_factorize(sparse, factor, &com);
   //
   // check properties of factor
   ok &= ( flag          == CHOLMOD_TRUE );  // return flag OK
   ok &= ( factor->n     == n );             // number of rows and coluns
   ok &= ( factor->minor == n );             // successful factorization
   ok &= ( factor->is_ll == CHOLMOD_FALSE ); // factorization is LDLT
   ok &= ( com.status == CHOLMOD_OK  );      // no problem with factorization
   //
   // information that defines the factor
   int*    L_perm = reinterpret_cast<int*>( factor->Perm );
   int*    L_nz   = reinterpret_cast<int*>( factor->nz );
   int*    L_p    = reinterpret_cast<int*>( factor->p );
   int*    L_i    = reinterpret_cast<int*>( factor->i );
   double* L_x    = reinterpret_cast<double*>( factor->x );
   //
   // check relation between Lz and Lp
   ok &= ( size_t(L_p[n]) == triplet->nnz );
   for(size_t j = 0; j < n; ++j)
      ok &= (L_nz[j] == L_p[j+1] - L_p[j] );
   //
   // P: permutation matrix corresponding to L_perm
   // multiplication by P on the right maps row i to row L_perm[i]
   std::vector<double> P(n*n, 0.0);
   for(size_t i = 0; i < n; ++i)
      P[L_perm[i] * n + i] = 1.0;

   //
   // L, D
   std::vector<double> L(n*n, 0.0), D(n*n, 0.0);
   for(size_t j = 0; j < n; ++j)
   {  for(size_t m = L_p[j]; m < size_t(L_p[j+1]); ++m)
      {  size_t i     = size_t( L_i[m] );
         if( i == j )
         {  D[i * n + i] = L_x[m];
            L[i * n + i] = 1.0;
         }
         else
            L[i * n + j] = L_x[m];
      }
   }
   //
   // LDLT = L * D * L^T
   std::vector<double> LD(n*n), LDLT(n*n);
   matrix_prod(n, false, L.data(), false, D.data(), LD.data());
   matrix_prod(n, false, LD.data(), true, L.data(), LDLT.data());
   //
   // PTAP = P^T * A * P
   std::vector<double> PTA(n*n), PTAP(n*n);
   matrix_prod(n, true, P.data(), false, A, PTA.data());
   matrix_prod(n, false, PTA.data(), false, P.data(), PTAP.data());
   //
   // check that P^T * A * P = L * D * L^T
# ifndef NDEBUG
   double eps99 = 99.0 * std::numeric_limits<double>::epsilon();
   for(size_t k = 0; k < n * n; ++k)
   {  double abs_err = std::fabs( PTAP[k] - LDLT[k] );
      assert( abs_err < eps99 );
   }
# endif
   //
   // free memory
   cholmod_free_triplet(&triplet,    &com);
   cholmod_free_sparse( &sparse,    &com);
   cholmod_free_factor( &factor,    &com);
   //
   // finish up
   cholmod_finish(&com);
   //
   // check count of malloc'ed - free'd
   ok &= com.malloc_count == 0;
   //
   return ok;
}
// END C++
