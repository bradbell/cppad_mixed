// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin ldlt_cholmod_update dev}
{xrst_spell
   factorize
   sym
}

Update Factorization Using new Matrix Values
############################################

Syntax
******
*ok* = *ldlt_obj* . ``update`` ( *H_rcv* )

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
This routine updates the :ref:`ldlt_cholmod-name` factorization
for new values in the square symmetric matrix.

ldlt_obj
********
This object has prototype

   ``CppAD::mixed::ldlt_cholmod`` *ldlt_obj*

In addition, it must have a previous call to
:ref:`ldlt_cholmod_init-name` .

H_rcv
*****
This argument contains new values for the
:ref:`sparse_mat_info@Notation@Sparse Matrix`
we are computing the LDLT factor of.
The :ref:`sparse_mat_info@Notation@Sparsity Pattern`
must be the same as in :ref:`ldlt_cholmod_init<ldlt_cholmod_init@H_rc>` .
Hence, in particular, it must be in
:ref:`column major<sparse_mat_info@Notation@Column Major Order>` order
and
:ref:`sparse_mat_info@Notation@Lower Triangular` .

sym_matrix\_
************
On input, the member variable

   ``cholmod_sparse`` * ``sym_matrix_``

has been :ref:`initialized<ldlt_cholmod_init@sym_matrix_>`
using the sparsity pattern as is in *H_rcv* .
Upon return, values in ``sym_matrix_`` have been set to the
corresponding values in the vector *H_rcv* . ``val`` () .

factor\_
********
On input, the member variable

   ``cholmod_factor`` * ``factor_``

has been :ref:`initialized<ldlt_cholmod_init@factor_>`
using the sparsity pattern for the Hessian.
Upon return, it contains the factorization

   ``cholmod_factorize`` ( ``sym_matrix_`` , ``factor_`` , & ``common_`` )

ok
**
If *ok* is true, the matrix was factored.
Otherwise, the matrix is singular.

Order of Operations
*******************
This *ldlt_obj* function must be called,
after the constructor and :ref:`init<ldlt_cholmod_init-name>`
and before any other member functions
(except :ref:`pattern<ldlt_cholmod_pattern-name>` ).

Example
*******
The file :ref:`ldlt_cholmod.cpp<ldlt_cholmod.cpp@update>` contains an
example and test that uses this function.

{xrst_end ldlt_cholmod_update}
*/
// ----------------------------------------------------------------------------


# include <cppad/mixed/include_cholmod.hpp>
# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cassert>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN_PROTOTYPE
bool ldlt_cholmod::update(const CppAD::mixed::d_sparse_rcv& H_rcv)
// END_PROTOTYPE
{  assert( init_done_ );
   update_called_ = true;
   //
   // set the values in sym_matrix_
   double* H_x  = reinterpret_cast<double*>( sym_matrix_->x);
# ifndef NDEBUG
   int*    H_p  = reinterpret_cast<int *>( sym_matrix_->p );
   int*    H_i  = reinterpret_cast<int *>( sym_matrix_->i );
# endif
   for(size_t k = 0; k < H_rcv.nnz(); k++)
   {  size_t ell = H_rc2cholmod_order_[k];
      H_x[ell]   = H_rcv.val()[k];
# ifndef NDEBUG
      assert( size_t( H_i[ell] ) == H_rcv.row()[k] );
      size_t j = H_rcv.col()[k];
      assert( size_t(H_p[j]) <= ell && ell < size_t(H_p[j+1]) );
# endif
   }
   // set factor_ to LDL^T factorization for this value of Hessian
# ifndef NDEBUG
   int flag =
# endif
   cholmod_factorize(sym_matrix_, factor_, &common_);
   //
   if( common_.status == CHOLMOD_NOT_POSDEF )
      return false;
   else
      assert( common_.status == CHOLMOD_OK );

   // check assumptions
   assert( flag           == CHOLMOD_TRUE );
   assert( factor_->n     == nrow_ );
   assert( factor_->minor <= nrow_ );
   assert( factor_->is_ll == CHOLMOD_FALSE );
   assert( factor_->xtype == CHOLMOD_REAL );
   // ---------------------------------------------------------------------
   if( factor_->minor < nrow_ )
      return false;
   return true;
}

} } // END_CPPAD_MIXED_NAMESPACE
