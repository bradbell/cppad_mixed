// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin ldlt_cholmod_inv dev}
{xrst_spell
   nrow
   sparseinv
}

Compute a Subset of the Inverse of Factored Matrix
##################################################

Under Construction
******************

Syntax
******
*ldlt_obj* . ``inv`` ( *row_in* , *col_in* , *val_out* )

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE
   // END_PROTOTYPE
}

Private
*******
The :ref:`ldlt_cholmod-name` class is an
:ref:`implementation detail<ldlt_cholmod@Private>` and not part of the
CppAD Mixed user API.

Purpose
*******
This function solves for a subset of the inverse of the
sparse symmetric matrix that has been factored.

ldlt_obj
********
This object has prototype

   ``const CppAD::mixed::ldlt_cholmod`` *ldlt_obj*

In addition, it must have a previous call to
:ref:`ldlt_cholmod_update-name` .

row_in
******
This vector contains the row indices for the components of the
inverse that we are computing.
It is assumed that these indices are the same for every call to
``inv`` .

col_in
******
This vector contains the column indices for the components of the
inverse that we are computing. It must have the same size as *row_in* .
It is assumed that these indices are the same for every call to
``inv`` .

val_out
*******
This matrix must have the same size as *row_in* .
The input values of its elements do not matter.
Upon return, it
contains the values for the components of the inverse that we are computing.
To be specific, for *k* = 0 , ..., *K* ``-1`` ,
*val_out* [ *k* ]
is *row_in* [ *k* ] , *col_in* [ *k* ] component of the inverse.

Method
******
This routine uses the
:ref:`install_unix@System Requirements@suitesparse` routine

   ``MATLAB_Tools/sparseinv/sparseinv.c``

to solve for the requested components of the inverse.

Member variables
****************
The first time this routine is called, the following member variable
as set to their permanent values:

sparseinv_p\_
=============
This member variable has the prototype

   ``CppAD::vector<int> sparseinv_p_``

Its initial size is zero.
After the first call to ``inv`` , it had size ``nrow_+1`` .
For *j* = 0, ..., ``nrow_`` ,
``sparseinv_p_`` [ *j* ] is where the row indices for the *j*-th
column start.

sparseinv_i\_
=============
This member variable has the prototype

   ``CppAD::vector<int> sparseinv_i_``

Its initial size is zero.
After the first call to ``inv`` , it had size

   ``sparseinv_p_`` [ ``nrow_`` ]

For each column index *j* , and for

   *k* = ``sparseinv_p_`` [ *j* ] , ..., ``sparseinv_p`` [ *j* +1] ``-1``

``sparseinv_i_`` [ *k* ] is the index of the next entry in the *j*-th
column of the sparse matrix.
Note that for each column, the corresponding row indices are in order; i.e.
sorted.
This sparse matrix includes all the indices in
*row_in* , *col_in* ,
all the indices in the factorization ``factor_`` ,
and the transpose of the indices in ``factor_`` .

sparseinv_order\_
=================
This member variable has the prototype

   ``CppAD::vector<size_t> sparseinv_order_``

Its initial size is zero.
After the first call to ``inv`` ,
it has the same size as *val_out* .
It is a mapping
from the indices in *row_in* , *col_in* and *val_out*
to the corresponding index in ``sparseinv_i`` .
Note that all the indices in
*row_in* and *col_in* also appear in the representation
``sparseinv_p_`` , ``sparseinv_i_`` .

{xrst_end ldlt_cholmod_inv}
------------------------------------------------------------------------------
*/
# include <cstddef>
# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cppad/mixed/sparseinv.hpp>
# include <cppad/utility/index_sort.hpp>

namespace CppAD { namespace mixed {

// BEGIN_PROTOTYPE
void ldlt_cholmod::inv(
   const CppAD::vector<size_t>& row_in    ,
   const CppAD::vector<size_t>& col_in    ,
   CppAD::vector<double>&       val_out   )
// END_PROTOTYPE
{  using CppAD::vector;
   //
   assert( update_called_ );
   assert( col_in.size() == row_in.size() );
   assert( val_out.size() == row_in.size() );
   //
   assert( factor_->n        == nrow_ );
   assert( factor_->minor    == nrow_ );
   assert( factor_->is_ll    == CHOLMOD_FALSE );
   assert( factor_->is_super == CHOLMOD_FALSE );
   assert( factor_->xtype    == CHOLMOD_REAL );
   //
   // information in the factor
   int*    L_perm  = reinterpret_cast<int*>( factor_->Perm );
   int*    L_p     = reinterpret_cast<int*>( factor_->p );
   int*    L_i     = reinterpret_cast<int*>( factor_->i );
   double* L_x     = reinterpret_cast<double*>( factor_->x );
   //
   // check order assumptions on factor
# ifndef NDEBUG
   int*    L_nz    = reinterpret_cast<int*>( factor_->nz );
   for(int j = 0; j < int( nrow_) ; j++)
   {  // diagonal entry first
      assert( L_i[ L_p[j] ] == j );
      assert( L_nz[j] == L_p[j+1] - L_p[j] );
      for(int k = L_p[j] + 1; k < L_p[j+1]; k++)
      {  // for each column, rows are sorted
         assert( L_i[k-1] < L_i[k] );
      }
   }
# endif
   //
   // number of rows and columns in the factor
   int n  = static_cast<int>(nrow_);
   //
   // extract the diagonal of the factor
   vector<double> d(nrow_);
   for(size_t j = 0; j < nrow_; j++)
   {  d[j] = L_x[ L_p[j] ];
      assert( d[j] != 0.0 );
   }
   //
   // number of requested indices
   size_t K_nnz  = row_in.size();
   // --------------------------------------------------------------------
   if( sparseinv_p_.size() == 0 )
   {  assert( sparseinv_i_.size() == 0 );
      assert( sparseinv_order_.size() == 0 );
      //
      // inverse permutation
      vector<size_t> p_inv(nrow_);
# ifndef NDEBUG
      for(size_t j = 0; j < nrow_; ++j)
         p_inv[j] = nrow_;
      for(size_t j = 0; j < nrow_; ++j)
         p_inv[ L_perm[j] ] = j;
      for(size_t j = 0; j < nrow_; ++j)
         assert( p_inv[j] != nrow_ );
# else
      for(size_t j = 0; j < nrow_; ++j)
         p_inv[ L_perm[j] ] = j;
# endif
      //
      // permuted version of indices
      vector<size_t> K_row(K_nnz), K_col(K_nnz);
      for(size_t k = 0; k < K_nnz; k++)
      {  K_row[k] = static_cast<size_t>( p_inv[ row_in[k] ] );
         K_col[k] = static_cast<size_t>( p_inv[ col_in[k] ] );
      }
      //
      // indices in K or L
      vector<size_t> KL_row( K_row ), KL_col( K_col );
      for(size_t j = 0; j < nrow_; j++)
      {  for(int m = L_p[j]; m < L_p[j+1]; ++m)
         {  size_t i = static_cast<size_t>( L_i[m] );
            // L
            KL_row.push_back(i);
            KL_col.push_back(j);
            if( i != j )
            {  KL_row.push_back(j);
               KL_col.push_back(i);
            }
         }
      }
      //
      // column major order
      size_t KL_nnz = KL_row.size();
      vector<size_t> key(KL_nnz), ind(KL_nnz);
      for(size_t k = 0; k < KL_nnz; k++)
         key[k] = KL_row[k] + KL_col[k] * nrow_;
      CppAD::index_sort(key, ind);
      //
      // put sparsity pattern in sparseinv_p_, sparseinv_i_
      // also fill in sparseinv_order_
      sparseinv_p_.resize(nrow_+1);
      sparseinv_order_.resize( K_nnz );
# ifndef NDEBUG
      for(size_t k = 0; k < K_nnz; ++k)
         sparseinv_order_[k] = KL_nnz; // invalid value
# endif
      sparseinv_p_[0] = 0;
      {  size_t ell = ind[0];
         size_t i   = KL_row[ell];
         size_t j   = KL_col[ell];
         //
         size_t i_ell = i;
         size_t j_ell = j;
         if( ell < K_nnz )
         {  i_ell = KL_row[ell];
            j_ell = KL_col[ell];
            sparseinv_order_[ell] = sparseinv_i_.size();
         }
         assert( i_ell == 0 && j_ell == 0 );
         //
         sparseinv_i_.push_back( int(i) );
         for(size_t k = 1; k < KL_nnz; ++k)
         {  ell = ind[k];
            if( ell < K_nnz )
            {  // same entry must not appear twice in row_in, col_in
               i_ell = KL_row[ell];
               j_ell = KL_col[ell];
               if( i_ell == i && j_ell == j )
               {  // this is previous entry pused to sparseinv_i_.
                  sparseinv_order_[ell] = sparseinv_i_.size() - 1;
               }
               else
               {  // this will be next entry pushed to sparseinv_i_
                  sparseinv_order_[ell] = sparseinv_i_.size();
               }
            }
            if( KL_col[ell] == j+1 )
            {  sparseinv_p_[j+1] = int( sparseinv_i_.size() );
               i = KL_row[ell];
               j = KL_col[ell];
               sparseinv_i_.push_back( int(i) );
            }
            else
            {  assert( KL_col[ell] == j );
               assert( i <= KL_row[ell] );
               if( i < KL_row[ell] )
               {  i = KL_row[ell];
                  sparseinv_i_.push_back( int(i) );
               }
            }
         }
         sparseinv_p_[nrow_] = int( sparseinv_i_.size() );
# ifndef NDEBUG
      for(size_t k = 0; k < K_nnz; ++k)
         assert( sparseinv_order_[k] < KL_nnz );
# endif
      }
   }
   // ------------------------------------------------------------------------
   // place where inverse is returned
   vector<double> Z_x( sparseinv_i_.size() );
   //
   // work space
   vector<double> z(nrow_);
   vector<int> Zdiagp(nrow_), Lmunch(nrow_);
   //
   // compute the subset of the inverse of the permuted matrix
   sparseinv(
      n,
      L_p,
      L_i,
      L_x,
      d.data(),
      L_p,
      L_i,
      L_x,
      sparseinv_p_.data(),
      sparseinv_i_.data(),
      Z_x.data(),
      z.data(),
      Zdiagp.data(),
      Lmunch.data()
   );
   //
   // return the values
   for(size_t k = 0; k < K_nnz; k++)
      val_out[k] = Z_x[ sparseinv_order_[k] ];
   //
   return;
}

} }
