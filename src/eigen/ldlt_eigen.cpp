// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-25 Bradley M. Bell
// ----------------------------------------------------------------------------

# include <Eigen/Sparse>
# include <cppad/mixed/ldlt_eigen.hpp>
# include <cppad/mixed/triple2eigen.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
-------------------------------------------------------------------------------
{xrst_begin ldlt_eigen_ctor dev}
{xrst_spell
  ptr
}

Eigen LDLT Constructor
######################

Syntax
******
``CppAD::ldlt_eigen<`` *Double* > *ldlt_obj* ( *n_row* )

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE_CTOR
   // END_PROTOTYPE_CTOR
}

Private
*******
The :ref:`ldlt_eigen-name` class is an
:ref:`implementation detail<ldlt_eigen@Private>` and not part of the
CppAD Mixed user API.

n_row\_
*******
The argument *n_row*
is the number of rows in the symmetric matrix we will compute factor.
The member variable ``n_row_`` is set to this value.

ptr\_
*****
This member variable points to a newly constructed
``eigen_ldlt`` object.

{xrst_end ldlt_eigen_ctor}
*/
// BEGIN_PROTOTYPE_CTOR
template <typename Double>
ldlt_eigen<Double>::ldlt_eigen(size_t n_row)
// END_PROTOTYPE_CTOR
:
n_row_(n_row),
init_done_(false),
update_called_(false)
{  ptr_ = new eigen_ldlt; }

// destructor
template <typename Double>
ldlt_eigen<Double>::~ldlt_eigen(void)
{  delete ptr_; }

/*
------------------------------------------------------------------------------
{xrst_begin ldlt_eigen_init dev}

Initialize LDLT Factor for a Specific Sparsity Pattern
######################################################

Syntax
******

   *ldlt_obj* . ``init`` ( *H_rc* )

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE_INIT
   // END_PROTOTYPE_INIT
}

Private
*******
The ``ldlt_eigen`` class is an
:ref:`implementation detail<ldlt_eigen@Private>` and not part of the
CppAD Mixed user API.

ldlt_obj
********
This object has prototype

   ``CppAD::ldlt_eigen<`` *Double* > *ldlt_obj*

H_rc
****
This argument is a
:ref:`sparse_mat_info@Notation@Sparsity Pattern` for the
matrix that we will compute the LDLT factor of.
It is in
:ref:`column major<sparse_mat_info@Notation@Column Major Order>` order
and
:ref:`sparse_mat_info@Notation@Lower Triangular` .

H_rc\_
======
This member variable is set to a copy of *H_rc* .

Order of Operations
*******************
This *ldlt_obj* function must be called once,
after the constructor and before any other member functions.

Example
*******
The file :ref:`ldlt_eigen.cpp<ldlt_eigen.cpp@init>` contains an
example and test that uses this function.

{xrst_end ldlt_eigen_init}
*/
// BEGIN_PROTOTYPE_INIT
template <typename Double>
void ldlt_eigen<Double>::init(const sparse_rc& H_rc)
// END_PROTOTYPE_INIT
{  assert( ! init_done_ );
   //
   // H_rc_
   H_rc_ = H_rc;
   //
   CppAD::vector<Double> not_used(0);
   eigen_sparse hessian_pattern;
   CppAD::mixed::triple2eigen(
         hessian_pattern  ,
         n_row_           ,
         n_row_           ,
         H_rc.row()       ,
         H_rc.col()       ,
         not_used
   );
   // analyze the pattern for an LDLT factorization of
   // f_{u,u}(theta, u)
   ptr_->analyzePattern(hessian_pattern);
   //
   init_done_ = true;
}
/*
-------------------------------------------------------------------------------
{xrst_begin ldlt_eigen_pattern dev}

Update Factorization Using new Matrix Values
############################################

Syntax
******
*H_rc* = *ldlt_obj* . ``pattern`` ()

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE_PATTERN
   // END_PROTOTYPE_PATTERN
}

Private
*******
The :ref:`ldlt_eigen-name` class is an
:ref:`implementation detail<ldlt_eigen@Private>` and not part of the
CppAD Mixed user API.

ldlt_obj
********
This object has prototype

   ``CppAD::mixed::ldlt_eigen`` *ldlt_obj*

In addition, it must have a previous call to
:ref:`ldlt_eigen_init-name` .

H_rc
****
The return value is a copy of the sparsity pattern
:ref:`ldlt_eigen_init@H_rc` in the corresponding call to
*ldlt_obj* . ``init`` ( *H_rc* ) .

Order of Operations
*******************
This *ldlt_obj* function must be called,
after the constructor and :ref:`init<ldlt_eigen_init-name>` .

Example
*******
The file :ref:`ldlt_eigen.cpp<ldlt_eigen.cpp@pattern>` contains an
example and test that uses this function.

{xrst_end ldlt_eigen_pattern}
*/
// BEGIN_PROTOTYPE_PATTERN
template <typename Double>
const sparse_rc& ldlt_eigen<Double>::pattern(void) const
// END_PROTOTYPE_PATTERN
{  assert( init_done_ );
   return H_rc_;
}

/*
------------------------------------------------------------------------------
{xrst_begin ldlt_eigen_update dev}
{xrst_spell
  ptr
}

Update Factorization Using new Matrix Values
############################################

Syntax
******
*ok* = *ldlt_obj* . ``update`` ( *H_rcv* )

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE_UPDATE
   // END_PROTOTYPE_UPDATE
}

Private
*******
The :ref:`ldlt_eigen-name` class is an
:ref:`implementation detail<ldlt_eigen@Private>` and not part of the
CppAD Mixed user API.

Purpose
*******
This routine updates the :ref:`ldlt_eigen-name` factorization
for new values in the square positive definite matrix.

ldlt_obj
********
This object has prototype

   ``CppAD::ldlt_eigen<`` *Double* > *ldlt_obj*

In addition, it must have a previous call to
:ref:`ldlt_eigen_init-name` .

H_rcv
*****
This argument contains new values for the
:ref:`sparse_mat_info@Notation@Sparse Matrix`
we are computing the LDLT factor of.
The :ref:`sparse_mat_info@Notation@Sparsity Pattern`
must be the same as in :ref:`ldlt_eigen_init<ldlt_eigen_init@H_rc>` .
Hence, in particular, it must be in
:ref:`column major<sparse_mat_info@Notation@Column Major Order>` order
and
:ref:`sparse_mat_info@Notation@Lower Triangular` .

ptr\_
*****
On input, the member variable *ptr_*
has been :ref:`initialized<ldlt_eigen_init-name>`
using the sparsity pattern for the Hessian.
Upon return, it contains the factorization

   ``ptr_->factorize`` ( *hessian* )

where *hessian* is an ``eigen_sparse``
representation of the Hessian with values.

ok
**
If the return value *ok* is true, the matrix was factored.
Otherwise, the matrix is singular.

Order of Operations
*******************
This *ldlt_obj* function must be called,
after the constructor and :ref:`init<ldlt_eigen_init-name>`
and before any other member functions.

Example
*******
The file :ref:`ldlt_eigen.cpp<ldlt_eigen.cpp@update>` contains an
example and test that uses this function.

{xrst_end ldlt_eigen_update}
*/
// BEGIN_PROTOTYPE_UPDATE
template <typename Double>
bool ldlt_eigen<Double>::update(
   const CppAD::sparse_rcv<s_vector, v_vector>& H_rcv
)
// END_PROTOTYPE_UPDATE
{  assert( init_done_ );
# ifndef NDEBUG
   assert( H_rcv.nnz() == H_rc_.nnz() );
   for(size_t k = 0; k < H_rc_.nnz(); ++k)
   {  assert( H_rcv.row()[k] == H_rc_.row()[k] );
      assert( H_rcv.col()[k] == H_rc_.col()[k] );
   }
# endif
   //
   eigen_sparse hessian;
   CppAD::mixed::triple2eigen(
      hessian     ,
      n_row_      ,
      n_row_      ,
      H_rcv.row() ,
      H_rcv.col() ,
      H_rcv.val()
   );
   // LDLT factorization of for specified values of the Hessian
   // f_{u,u}(theta, u)
   ptr_->factorize(hessian);
   //
   if( ptr_->info() != Eigen::Success )
      return false;
   //
   update_called_ = true;
   return true;
}
/*
-------------------------------------------------------------------------------
{xrst_begin ldlt_eigen_split dev}

Split Out a Factorization
#########################

Syntax
******
*ldlt_obj* . ``split`` ( *L* , *D* , *P* )

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE_SPLIT
   // END_PROTOTYPE_SPLIT
}

Private
*******
The :ref:`ldlt_eigen-name` class is an
:ref:`implementation detail<ldlt_eigen@Private>` and not part of the
CppAD Mixed user API.

Purpose
*******
Extract the components of the :ref:`ldlt_eigen@Factorization`

.. math::

   L D L^\R{T} = P H P^{T}

ldlt_obj
********
This object has prototype

   ``const CppAD::ldlt_eigen<`` *Double* > *ldlt_obj*

In addition, it must have a previous call to
:ref:`ldlt_eigen_update-name` .

L
*
is a lower triangular matrix with ones on the diagonal.

D
*
is a diagonal matrix.

P
*
is a permutation matrix.

Example
*******
The file :ref:`ldlt_eigen.cpp<ldlt_eigen.cpp@split>` contains an
example and test that uses this function.

{xrst_end ldlt_eigen_split}
*/

// BEGIN_PROTOTYPE_SPLIT
template <typename Double>
void ldlt_eigen<Double>::split(
   Eigen::SparseMatrix<Double, Eigen::ColMajor>& L ,
   Eigen::Matrix<Double, Eigen::Dynamic, 1>&     D ,
   Eigen::PermutationMatrix<Eigen::Dynamic>&     P ) const
// END_PROTOTYPE_SPLIT
{  assert( update_called_ );
   //
   L = ptr_->matrixL();
   D = ptr_->vectorD();
   P = ptr_->permutationP();
}
/*
------------------------------------------------------------------------------
{xrst_begin ldlt_eigen_rcond dev}

Reciprocal of Condition Number for D
####################################

Syntax
******
*rcond* = *ldlt_obj* . ``rcond`` ()

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE_RCOND
   // END_PROTOTYPE_RCOND
}

Private
*******
The ``ldlt_eigen`` class is an
:ref:`implementation detail<ldlt_eigen@Private>` and not part of the
CppAD Mixed user API.

ldlt_obj
********
The object *ldlt_obj*
must have a previous call to :ref:`ldlt_eigen_update-name` .


rcond
*****
The return value *rcond*
is the reciprocal of the condition number for the diagonal matrix *D*
in the factorization.
In other words, it is the minimum absolute entry in *D* divided
by the maximum absolute entry in *D* .
If the matrix *D* is singular, or any entry in *D* nan or infinite,
*rcond* is zero.

Example
*******
The file :ref:`ldlt_eigen.cpp<ldlt_eigen.cpp@rcond>` contains an
example and test that uses this function.

{xrst_end ldlt_eigen_rcond}
*/
// BEGIN_PROTOTYPE_RCOND
template <typename Double>
Double ldlt_eigen<Double>::rcond(void) const
// END_PROTOTYPE_RCOND
{  assert( update_called_ );

   // diag
   Eigen::Matrix<Double, Eigen::Dynamic, 1> diag = ptr_->vectorD();
   assert( diag.size() == int(n_row_) );
   //
   // max_abs, min_abs
   Double max_abs = 0.0;
   Double min_abs = CppAD::numeric_limits<Double>::infinity();
   for(size_t j = 0; j < n_row_; j++)
   {  Double abs = fabs( diag[j] );
      if( isnan( abs ) )
         abs = 0.0;
      max_abs    = std::max( abs, max_abs);
      min_abs    = std::min( abs, min_abs);
   }
   //
   // rcond
   if( min_abs == 0.0 )
      return 0.0;
   if( min_abs == CppAD::numeric_limits<Double>::infinity() )
      return 0.0;
   if( max_abs == CppAD::numeric_limits<Double>::infinity() )
      return 0.0;
   return min_abs / max_abs;
}
/*
------------------------------------------------------------------------------
{xrst_begin ldlt_eigen_logdet dev}

Compute Log Determinant for Current LDLT Factor
###############################################

Syntax
******
*logdet* = *ldlt_obj* . ``logdet`` ( *negative* )

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE_LOGDET
   // END_PROTOTYPE_LOGDET
}

Private
*******
The ``ldlt_eigen`` class is an
:ref:`implementation detail<ldlt_eigen@Private>` and not part of the
CppAD Mixed user API.

ldlt_obj
********
This object has prototype

   ``const CppAD::ldlt_eigen<`` *Double* > *ldlt_obj*

In addition, it must have a previous call to
:ref:`ldlt_eigen_update-name` .

negative
********
The input value of *negative* does no matter,
upon return it is the number of elements of
:ref:`ldlt_eigen@Factorization@D`
that are less than zero.

logdet
******
This return value *logdet*
is the log of the absolute value of the determinant corresponding
to the previous call to :ref:`ldlt_eigen_update-name` .
If the matrix is singular, *logdet* is
minus infinity.

Example
*******
The file :ref:`ldlt_eigen.cpp<ldlt_eigen.cpp@logdet>` contains an
example and test that uses this function.

{xrst_end ldlt_eigen_logdet}
*/
// BEGIN_PROTOTYPE_LOGDET
template <typename Double>
Double ldlt_eigen<Double>::logdet(size_t& negative) const
// END_PROTOTYPE_LOGDET
{  assert( update_called_ );

   // compute the logdet( f_{u,u}(theta, u )
   Eigen::Matrix<Double, Eigen::Dynamic, 1> diag = ptr_->vectorD();
   assert( diag.size() == int(n_row_) );
   negative        = 0;
   bool   has_zero = false;
   Double logdet   = 0.0;
   for(size_t j = 0; j < n_row_; j++)
   {  has_zero |= diag(j) == 0.0;
      if( diag(j) < 0.0 )
         negative++;
      logdet += log( CppAD::fabs( diag(j) ) );
   }
   if( has_zero )
      return - std::numeric_limits<Double>::infinity();
   //
   return logdet;
}
/*
-----------------------------------------------------------------------------
{xrst_begin ldlt_eigen_solve_H dev}

Solve Linear Equations Using Stored Factor
##########################################

Syntax
******
*ldlt_obj* . ``solve_H`` ( *row* , *val_in* , *val_out* )

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE_SOLVE_H
   // END_PROTOTYPE_SOLVE_H
}

Private
*******
The :ref:`ldlt_eigen-name` class is an
:ref:`implementation detail<ldlt_eigen@Private>` and not part of the
CppAD Mixed user API.

Purpose
*******
This function solves the linear equation
:math:`H x = b` where :math:`H` is the positive definite matrix
corresponding to the previous :ref:`update<ldlt_eigen_update-name>` ,
:math:`b` is a known column vector,
and :math:`x` is unknown.

ldlt_obj
********
This object has prototype

   ``const CppAD::ldlt_eigen<`` *Double* > *ldlt_obj*

In addition, it must have a previous call to
:ref:`ldlt_eigen_update-name` .

row
***
This argument
contains all of the rows of column vector :math:`b` that are
non-zero and the rows of the column vector *x*
that are desired.
These values are in strictly increasing order; i.e.,

   *row* [ *k* ] < *row* [ *k* +1]

It follows that *row* . ``size`` () is less than or equal
:ref:`ldlt_eigen_ctor@n_row_` .

val_in
******
This argument has the same size as *row* .
It specifies the values in the column vector :math:`b`
for each of the corresponding rows; i.e.,
for *k* = 0 , ..., *row* . ``size`` () ``-1`` ,

   *b* [ *row* [ *k* ] ] = *val_in* [ *k* ]

.

val_out
*******
This argument has the same size as *row* .
On input, the value of its elements do not matter.
Upon return, it contains the values in the column vector :math:`b`
for each of the corresponding rows; i.e.,
for *k* = 0 , ..., *row* . ``size`` () ``-1`` ,

   *x* [ *row* [ *k* ] ] = *val_out* [ *k* ]

.

Example
*******
The file :ref:`ldlt_eigen.cpp<ldlt_eigen.cpp@solve_H>` contains an
example and test that uses this function.

{xrst_end ldlt_eigen_solve_H}
*/
// BEGIN_PROTOTYPE_SOLVE_H
template <typename Double>
void ldlt_eigen<Double>::solve_H(
   const s_vector&              row     ,
   const CppAD::vector<Double>& val_in  ,
   CppAD::vector<Double>&       val_out ) const
// END_PROTOTYPE_SOLVE_H
{  assert( update_called_ );
   assert( row.size() == val_in.size() );
   assert( row.size() == val_out.size() );
   //
   // eigen uses dense vectors during a sparse solve so do the same
   // (in fact, sparse vectors here fails some test as if short circuiting)
   eigen_vector b = eigen_vector::Zero(n_row_);
   for(size_t k = 0; k < row.size(); k++)
   {  assert( row[k] < n_row_ );
      b[ row[k]] = val_in[k];
   }
   //
   // val_out
   eigen_vector x = ptr_->solve(b);
   for(size_t k = 0; k < row.size(); k++)
      val_out[k] = x[ row[k] ];
}
/*
-------------------------------------------------------------------------------
{xrst_begin ldlt_eigen_sim_cov dev}
{xrst_spell
  ll
}

Simulations with Covariance Corresponding to Factored Matrix
############################################################

Syntax
******
*ok* = *ldlt_obj* . ``sim_cov`` ( *w* , *v* )

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE_SIM_COV
   // END_PROTOTYPE_SIM_COV
}

Private
*******
The :ref:`ldlt_eigen-name` class is an
:ref:`implementation detail<ldlt_eigen@Private>` and not part of the
CppAD Mixed user API.

Purpose
*******
This function simulates a normal random vector with mean zero
and covariance :math:`H^{-1}` where

.. math::

   L D L^\R{T} = P H P^\R{T}

is the current factorization; see
:ref:`ldlt_eigen@Factorization@H` ,
:ref:`ldlt_eigen@Factorization@L` ,
:ref:`ldlt_eigen@Factorization@D` , and
:ref:`ldlt_eigen@Factorization@P` .

ldlt_obj
********
This object has prototype

   ``const CppAD::ldlt_eigen<`` *Double* > *ldlt_obj*

In addition, it must have a previous call to
:ref:`ldlt_eigen_update-name` .

w
*
This argument's
size is equal to the number of rows in :math:`H`.

v
*
This argument's
size is equal to the number of rows in :math:`H`.
The input value of its elements does not matter.
Upon return

.. math::

   v = P^\R{T} L^{-\R{T}} \tilde{D}^{-1/2} w

If :math:`w` is mean zero, variance identity white noise,
:math:`w \sim \B{N} ( 0 , I )`,
then :math:`v` will be mean zero and variance :math:`H^{-1}`,
:math:`v \sim \B{N} ( 0 , H^{-1} )`; see
:ref:`theory@Sparse Observed Information` .

Positive Definite
*****************
In the formula for :math:`v` above,
the matrix :math:`\tilde{D}` is a positive version of :math:`D`.
To be specific,

.. math::

   \tilde{D}_{i,i} = \left\{ \begin{array}{ll}
      D_{i,i} & \R{if} \; D_{i,i} \geq  \varepsilon^2 \; \max(D) \\
      \varepsilon^2 \; \max(D) & \R{otherwise}
   \end{array} \right.

where :math:`\varepsilon`
``std::numeric_limits<Double>::epsilon()`` ,
and :math:`\max(D)` is the largest element in :math:`D`.

ok
**
The return value has prototype

   ``bool`` *ok*

If :math:`\max(D) > 0`, this routine terminates with *ok*
equal to true.
Otherwise it is false and the output values in *v*
are the same as their input values.

Example
*******
The file :ref:`ldlt_eigen.cpp<ldlt_eigen.cpp@sim_cov>` contains an
example and test that uses this function.

{xrst_end ldlt_eigen_sim_cov}
*/
// BEGIN_PROTOTYPE_SIM_COV
template <typename Double>
bool ldlt_eigen<Double>::sim_cov(
   const CppAD::vector<Double>& w  ,
   CppAD::vector<Double>&       v  ) const
// END_PROTOTYPE_SIM_COV
{  assert( update_called_ );
   typedef Eigen::Matrix<Double, Eigen::Dynamic, 1> column_vector;
   //
   // set b = w
   column_vector b(n_row_);
   for(size_t i = 0; i < n_row_; i++)
      b[i] = w[i];
   //
   // diagonal
   column_vector diag = ptr_->vectorD();
   Double max_D = 0.0;
   for(size_t i = 0; i < n_row_; i++)
      max_D = std::max(max_D, diag[i] );
   if( max_D <= 0.0 )
      return false;
   //
   // set b = D^{-1/2} w
   Double eps = std::numeric_limits<Double>::epsilon();
   eps        = eps * eps * max_D;
   for(size_t i = 0; i < n_row_; i++)
   {  Double di = std::max(diag[i], eps);
      b[i] = b[i] / CppAD::sqrt( di );
   }
   //
   // set b = L^{-T} * D^{-1/2} w
   b = ptr_->matrixU().solve(b);
   //
   // set b = P^T L^{-T} * D^{-1/2} w
   b = ptr_->permutationP().transpose() * b;
   //
   // return v
   for(size_t i = 0 ; i < n_row_; i++)
      v[i] = b[i];
   //
   return true;
}

/*
-------------------------------------------------------------------------------
{xrst_begin ldlt_eigen_inv dev}

Compute a Subset of the Inverse of Factored Matrix
##################################################

Syntax
******
*ldlt_obj* . ``inv`` ( *row_in* , *col_in* , *val_out* )

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE_INV
   // END_PROTOTYPE_INV
}

Private
*******
The :ref:`ldlt_eigen-name` class is an
:ref:`implementation detail<ldlt_eigen@Private>` and not part of the
CppAD Mixed user API.

Purpose
*******
This function solves for a subset of the inverse of the
sparse symmetric matrix that has been factored.

ldlt_obj
********
This object has prototype

   ``const CppAD::ldlt_eigen<`` *Double* > *ldlt_obj*

In addition, it must have a previous call to
:ref:`ldlt_eigen_update-name` .

row_in
******
This vector contains the row indices for the components of the
inverse that we are computing.

col_in
******
This vector contains the column indices for the components of the
inverse that we are computing. It must have the same size as *row_in* .

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
This routine uses :ref:`ldlt_eigen_solve_H-name` to solve for the
requested components of the inverse one column at a time.

{xrst_end ldlt_eigen_inv}
*/

// BEGIN_PROTOTYPE_INV
template <typename Double>
void ldlt_eigen<Double>::inv(
   const s_vector&              row_in    ,
   const s_vector&              col_in    ,
   CppAD::vector<Double>&       val_out   ) const
// END_PROTOTYPE_INV
{  assert( update_called_ );
   using CppAD::vector;
   //
   // starting index
   size_t K  = row_in.size();
   size_t k  = 0;
   size_t row = n_row_;
   size_t col = n_row_;
   if( k < K )
   {  row = row_in[k];
      col = col_in[k];
      assert( row < n_row_ );
      assert( col < n_row_ );
   }
   //
   s_vector row_solve;
   CppAD::vector<Double> rhs_solve, val_solve;
   for(size_t j = 0; j < n_row_; j++)
   {  // vectors for this column
      row_solve.resize(0);
      rhs_solve.resize(0);

      // only need for rows where f_{u,u} (theta_u) is possibly not zero
      size_t k_start        = K;
      bool   found_diagonal = false;
      while( col <= j )
      {  if( col == j )
         {  // this row of f_{u,u} (theta, u) is possibly non-zero
            row_solve.push_back( row );
            // rhs_solve needs to be j-th column of identity matrix
            // so solution is j-th column of inverse
            if( row == j )
            {  rhs_solve.push_back(1.0);
               found_diagonal = true;
            }
            else
               rhs_solve.push_back(0.0);
            //
            // index in row_in and col_in where j-th column starts
            if( k_start == K )
               k_start = k;
         }
         k++;
         if( k < K )
         {  row = row_in[k];
            col = col_in[k];
            assert( row < n_row_ );
            assert( col < n_row_ );
         }
         else
            row = col = n_row_;
      }
      assert( col > j );
      //
      // Cannot compute cholesky factor if f_{u,u} (theta, u) is zero
      if( ! found_diagonal )
      {  row_solve.push_back(j);
         rhs_solve.push_back(1.0);
      }
      // if k_start == K, we do not need any components of the inverse
      if( k_start < K )
      {  val_solve.resize( row_solve.size() );
         //
         solve_H(row_solve, rhs_solve, val_solve);
         //
         size_t nr = row_solve.size();
         if( ! found_diagonal )
            nr--;
         for(size_t ell = 0; ell < nr; ell++)
            val_out[k_start + ell] = val_solve[ell];
      }
   }
   return;
}
/*
-------------------------------------------------------------------------------
{xrst_begin ldlt_eigen_solve_LDLT dev}

Solve Linear Equations Corresponding to L, D, and P Factors
###########################################################

Syntax
******
*x* = ``ldlt_eigen`` < *Double* >:: ``solve_LDLT`` ( *L* , *D* , *P* , *b* )

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE_SOLVE_LDLT
   // END_PROTOTYPE_SOLVE_LDLT
}

Private
*******
The :ref:`ldlt_eigen-name` class is an
:ref:`implementation detail<ldlt_eigen@Private>` and not part of the
CppAD Mixed user API.

Purpose
*******
This function solves the linear equation
:math:`H x = b` where the positive definite matrix
:math:`H = P^T L D L^T P`,
:math:`b` is a known column vector,
and :math:`x` is unknown.

L
*
is a lower triangular matrix with ones on the diagonal.

D
*
is a diagonal matrix.

P
*
is a permutation matrix.

b
*
is the right hand side column vector in the equation.

x
*
is the column vector that solves the equation.

Example
*******
The file :ref:`ldlt_eigen.cpp<ldlt_eigen.cpp@solve_LDLT>` contains an
example and test that uses this function.

{xrst_end ldlt_eigen_solve_LDLT}
*/
// BEGIN_PROTOTYPE_SOLVE_LDLT
template <typename Double>
Eigen::Matrix<Double, Eigen::Dynamic, 1>
ldlt_eigen<Double>::solve_LDLT(
   // Lower Triangular matrix
   const Eigen::SparseMatrix<Double, Eigen::ColMajor>& L  ,
   // Diagonal matrix
   const eigen_vector&                                 D  ,
   // Permutation matrix
   const Eigen::PermutationMatrix<Eigen::Dynamic>&     P  ,
   // right hand side of equation
   const eigen_vector&                                 b  )
// END_PROTOTYPE_SOLVE_LDLT
{  using Eigen::Lower;
   using Eigen::Upper;
   //
   size_t n = size_t( b.size() );
   //
   // P * b
   eigen_vector result = P * b;
   //
   // L^-1 * P * b
   result = L. template triangularView<Lower>().solve(result);
   //
   // D^-1 * L^-1 * P * b
   for(size_t j = 0; j < n; ++j)
      result[j] = result[j] / D[j];
   //
   // L^-T * D^-1 * L^-1 * P * b
   result = L.transpose(). template triangularView<Upper>().solve(result);
   //
   // P^T * L^-T * D^-1 * L^-1 * P * b
   result = P.transpose() * result;
   //
   return result;
}

} } // END_CPPAD_MIXED_NAMESPACE

// Explicit instantiation of ldlt_eigen
template class CppAD::mixed::ldlt_eigen<double>;
template class CppAD::mixed::ldlt_eigen< CppAD::mixed::a1_double >;
