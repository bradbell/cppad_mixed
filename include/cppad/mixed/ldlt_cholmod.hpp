// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-25 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_LDLT_CHOLMOD_HPP
# define CPPAD_MIXED_LDLT_CHOLMOD_HPP

/*
{xrst_begin ldlt_cholmod dev}
{xrst_spell
  preprocessor
}

A Cholmod Cholesky Factor Class
###############################

See Also
********
:ref:`ldlt_eigen-name`

Private
*******
This class is an implementation detail and not part of the
CppAD Mixed user API.

Purpose
*******
This class has utilities that work with a ``cholmod`` LDLT factor,
``cholmod`` is part of the
:ref:`install_unix@System Requirements@suitesparse` package.

Factorization
*************
The factorization is

.. math::

   L D L^\R{T} = P H P^{T}

where

H
=
is the matrix corresponding the current
:ref:`update<ldlt_cholmod_update-name>` .

L
=
is a lower triangular matrix with ones on the diagonal,

D
=
is a diagonal matrix.

P
=
is a permutation matrix.

Example
*******
The file :ref:`ldlt_cholmod.cpp-name` contains an example and test
using the operations in this class.

Preprocessor Symbols
********************
The following extra ``CHOLMOD_`` * symbols are defined
{xrst_literal
   // BEGIN SYMBOLS
   // END SYMBOLS
}

Contents
********
{xrst_toc_table
   src/cholmod/constructor.cpp
   src/cholmod/init.cpp
   src/cholmod/pattern.cpp
   src/cholmod/update.cpp
   src/cholmod/rcond.cpp
   src/cholmod/logdet.cpp
   src/cholmod/solve_H.cpp
   src/cholmod/sim_cov.cpp
   src/cholmod/inv.cpp
   example/private/ldlt_cholmod.cpp
   example/private/cholmod_factor.cpp
   example/private/cholmod_solve.cpp
   example/private/cholmod_solve2_a.cpp
   example/private/cholmod_solve2_sim.cpp
   example/private/sparseinv.cpp
}

{xrst_end ldlt_cholmod}
------------------------------------------------------------------------------
*/

# include <cppad/mixed/include_cholmod.hpp>
# include <cppad/mixed/typedef.hpp>

// BEGIN SYMBOLS
# define CHOLMOD_TRUE                  1
# define CHOLMOD_FALSE                 0
# define CHOLMOD_STYPE_NOT_SYMMETRIC   0
# define CHOLMOD_STYPE_LOWER_TRIANGLE -1
// END SYMBOLS


namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

class ldlt_cholmod {
private:
   //
   const size_t    nrow_;           // number of rows (and columns) in H
   bool            init_done_;     // has init been called
   bool            update_called_; // has update been called
   sparse_rc       H_rc_;          // sparsity pattern for H
   //
   CppAD::vector<size_t> key_;    // temporary sorting keys
   CppAD::vector<size_t> index_;  // temporary sorting indices
   //
   // mapping from sparse_mat_info to cholmod_sparse order
   CppAD::vector<size_t> H_rc2cholmod_order_;
   //
   cholmod_common    common_;
   cholmod_sparse*   sym_matrix_; // The symmetric matrix we are factoring
   cholmod_factor*   factor_;     // lower triangular LDL' factor
   //
   // information used by inv and computed during first call to inv
   CppAD::vector<size_t> sparseinv_order_;
   CppAD::vector<int>    sparseinv_p_, sparseinv_i_;
   //
   //
   // work space for solving equations
   cholmod_dense*    rhs_;        // right hand side of equation
   cholmod_sparse*   rhs_set_;    // sparsity pattern for rhs
   cholmod_dense*    sol_;        // solution of equaiton
   cholmod_sparse*   sol_set_;    // sparsity pattern for solutions
   cholmod_dense*    work_one_;   // first work space for solving equations
   cholmod_dense*    work_two_;   // second work space for solving equations
public:
   // ----------------------------------------------------------------------
   // non-const functions
   //
   // constructor
   ldlt_cholmod(size_t nrow);
   // destructor
   ~ldlt_cholmod(void);
   // initialize
   void init(const sparse_rc& hes_rc);
   // factorize
   bool update(const d_sparse_rcv& hes_rcv);
   // ----------------------------------------------------------------------
   // const functions
   //
   // pattern
   const sparse_rc& pattern(void) const;
   //
   // rcond
   // reciprocal of condition number of D
   double rcond(void) const;
   //
   // logdet
   // log determinant
   double logdet(size_t& negative) const;
   //
   // compute a subset of the inverse
   void inv(
      const CppAD::vector<size_t>& row      ,
      const CppAD::vector<size_t>& col      ,
      CppAD::vector<double>&       val
   );
   // solve linear equations
   void solve_H(
      const CppAD::vector<size_t>& row      ,
      const CppAD::vector<double>& val_in   ,
      CppAD::vector<double>&       val_out
   );
   // simualte covariance
   bool sim_cov(
      const CppAD::vector<double>& w ,
      CppAD::vector<double>&       v
   );
};


} } // END_CPPAD_MIXED_NAMESPACE


# endif
