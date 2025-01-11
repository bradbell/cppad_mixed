// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin ldlt_cholmod_init dev}
{xrst_spell
  nrow
  rhs
  sym
}

Initialize Factor for a Specific Sparsity Pattern
#################################################

Syntax
******

   *ldlt_obj* . ``init`` ( *H_rc* )

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

ldlt_obj
********
This object has prototype

   ``CppAD::mixed::ldlt_cholmod`` *ldlt_obj*

H_rc
****
This argument is a
:ref:`sparse_mat_info@Notation@Sparsity Pattern` for the
square matrices with
:ref:`ldlt_cholmod_ctor@nrow_` rows that we will compute the
LDLT factor of.
It is in
:ref:`column major<sparse_mat_info@Notation@Column Major Order>` order
and
:ref:`sparse_mat_info@Notation@Lower Triangular` .

H_rc\_
======
This member variable is set to a copy of *H_rc* .

Assumptions
***********
All of the ``cholmod`` private pointers
are null when this routine is called.

sym_matrix\_
************
Upon return,
this is set to a packed, real, sorted, lower triangular, sparse matrix
with the pattern specified by *H_rc* and
the value ``nan`` for each possibly non-zero value.

factor\_
********
Upon return,
this is the result of

   ``factor_`` = ``cholmod_analyze`` ( ``sym_matrix_`` , & ``common_`` )

rhs\_
*****
Upon return,
this is set to dense column vector of ``nrow_`` zeros as follows:

   ``rhs_`` = ``cholmod_zeros`` ( ``nrow_`` , 1, ``CHOLMOD_REAL`` , & ``common_`` )

rhs_set\_
*********
Upon return,
this is a sparse column vector with ``nrow_`` rows.
It is packed, sorted, not symmetric, and just a sparsity pattern.
(Note that sorted refers the column order and does not matter for vectors.)
There are ``nrow_``
possibly non-zero elements, but which elements
are non-zero is not specified
(expected to be set for each :ref:`ldlt_cholmod_solve_H-name` usage).

Order of Operations
*******************
This *ldlt_obj* function must be called once,
after the constructor and before any other member functions.

Example
*******
The file :ref:`ldlt_cholmod.cpp<ldlt_cholmod.cpp@init>` contains an
example and test that uses this function.

{xrst_end ldlt_cholmod_init}
*/

# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cppad/utility/index_sort.hpp>
# include <cassert>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN_PROTOTYPE
void ldlt_cholmod::init(const CppAD::mixed::sparse_rc& H_rc)
// END_PROTOTYPE
{  assert( ! init_done_ );
   //
   // H_rc_
   H_rc_ = H_rc;
   //
   assert(sym_matrix_ == CPPAD_NULL );
   assert(factor_     == CPPAD_NULL );
   assert(rhs_        == CPPAD_NULL );
   assert(rhs_set_    == CPPAD_NULL );
   assert(sol_        == CPPAD_NULL );
   assert(sol_set_    == CPPAD_NULL );
   assert(work_one_   == CPPAD_NULL );
   assert(work_two_   == CPPAD_NULL );

   // set triplet corresponding to sparsity pattern and reserve spase
   // for unspecified values (that end up in sym_matrix_.
   size_t ncol  = nrow_;
   size_t nzmax = H_rc.nnz();
   int    stype = CHOLMOD_STYPE_LOWER_TRIANGLE;
   int    xtype = CHOLMOD_REAL;
   cholmod_triplet* triplet = cholmod_allocate_triplet(
      nrow_, ncol, nzmax, stype, xtype, &common_
   );
   int* T_i    = (int *) triplet->i;
   int* T_j    = (int *) triplet->j;
   double* T_x = (double *) triplet->x;
   double nan  = std::numeric_limits<double>::quiet_NaN();
   for(size_t k = 0; k < nzmax; k++)
   {  assert( H_rc.row()[k] < nrow_ );
      assert( H_rc.col()[k] < nrow_ );
      //
      T_i[k] = static_cast<int>( H_rc.row()[k] );
      T_j[k] = static_cast<int>( H_rc.col()[k] );
      //
      T_x[k] = nan;
   }
   triplet->nnz = nzmax;

   // convert triplet to a sparse matrix
   nzmax = 0; // just need max corresponding to the triplet
   sym_matrix_ = cholmod_triplet_to_sparse(triplet, nzmax, &common_);

   // check assumptions
   assert( sym_matrix_->nrow   == nrow_ );
   assert( sym_matrix_->ncol   == nrow_ );
   assert( sym_matrix_->nzmax  == H_rc.nnz() );
   assert( sym_matrix_->stype  == CHOLMOD_STYPE_LOWER_TRIANGLE );
   assert( sym_matrix_->itype  == CHOLMOD_INT );
   assert( sym_matrix_->xtype  == CHOLMOD_REAL );
   assert( sym_matrix_->dtype  == CHOLMOD_DOUBLE );
   assert( sym_matrix_->sorted == CHOLMOD_TRUE  );
   assert( sym_matrix_->packed == CHOLMOD_TRUE  );


   // analyze the sparsity pattern for LDLT factorization of
   factor_ = cholmod_analyze(sym_matrix_, &common_);

   // check assumptions
   assert( factor_->n     == nrow_ );
   assert( factor_->minor == nrow_ );
   assert( factor_->is_ll == CHOLMOD_FALSE );

   // set rhs_ to column vector of zeros
   rhs_ = cholmod_zeros(nrow_, 1, CHOLMOD_REAL, &common_);

   // set rhs_set_ to be a sparsity pattern, but do not specify the
   // actual non-zero entries.
   ncol       = 1;
   nzmax      = nrow_;
   int sorted = CHOLMOD_TRUE;
   int packed = CHOLMOD_TRUE;
   stype      = CHOLMOD_STYPE_NOT_SYMMETRIC;
   xtype      = CHOLMOD_PATTERN;
   rhs_set_   = cholmod_allocate_sparse(
      nrow_,
      ncol,
      nzmax,
      sorted,
      packed,
      stype,
      xtype,
      &common_
   );
   // check
   assert( rhs_set_->nrow   == nrow_ );
   assert( rhs_set_->ncol   == 1     );
   assert( rhs_set_->nzmax  == nrow_ );
   assert( rhs_set_->stype  == CHOLMOD_STYPE_NOT_SYMMETRIC );
   assert( rhs_set_->itype  == CHOLMOD_INT );
   assert( rhs_set_->sorted == CHOLMOD_TRUE  );
   assert( rhs_set_->packed == CHOLMOD_TRUE  );

   // done with triplet
   cholmod_free_triplet(&triplet, &common_ );
   //
   // sym_matrix_ is in column major order
   H_rc2cholmod_order_ = H_rc.col_major();
   //
# ifndef NDEBUG
   int*    H_p  = reinterpret_cast<int *>( sym_matrix_->p );
   int*    H_i  = reinterpret_cast<int *>( sym_matrix_->i );
   size_t  nnz  = H_rc.nnz();
   CppAD::vector<size_t> cholmod2H_rc_order(nnz);
   for(size_t k = 0; k < nnz; k++)
   {  size_t ell = H_rc2cholmod_order_[k];
      cholmod2H_rc_order[ell] = k;
   }
   size_t ell = 0;
   for(size_t c = 0; c < nrow_; c++)
   {  for(int m = H_p[c]; m < H_p[c+1]; ++m)
      {  size_t k  = cholmod2H_rc_order[ell];
         assert( H_rc.row()[k]  == size_t( H_i[m] ) );
         assert( H_rc.col()[k]  == c );
         ell++;
      }
   }
# endif
   init_done_ = true;
}

} } // END_CPPAD_MIXED_NAMESPACE
