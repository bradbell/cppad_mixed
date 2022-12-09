// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_LDLT_EIGEN_HPP
# define CPPAD_MIXED_LDLT_EIGEN_HPP

/*
{xrst_begin ldlt_eigen dev}

An Eigen LDLT Factor Class
##########################

See Also
********
:ref:`ldlt_cholmod-name`

Private
*******
This class is an implementation detail and not part of the
CppAD Mixed user API.

Factorization
*************
The factorization is

.. math::

   L D L^\R{T} = P H P^{T}

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

Double
******
This is the type of the elements in the matrices and is either
 double29628 or  a1_double29628.

Example
*******
The file :ref:`ldlt_eigen.cpp-name` contains an example and test
using the operations in this class.

Contents
********
{xrst_toc_table
   src/eigen/ldlt_eigen.cpp
   example/private/ldlt_eigen.cpp
}

{xrst_end ldlt_eigen}
------------------------------------------------------------------------------
*/

# include <Eigen/Sparse>
# include <cppad/cppad.hpp>
# include <cppad/mixed/typedef.hpp>


namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

template <typename Double>
class ldlt_eigen {
private:
   typedef CppAD::vector<Double>                             v_vector;
   typedef Eigen::SparseMatrix<Double, Eigen::ColMajor>      eigen_sparse;
   typedef Eigen::Matrix<Double, Eigen::Dynamic, 1>          eigen_vector;
   typedef Eigen::PermutationMatrix<Eigen::Dynamic>          eigen_perm;
   typedef Eigen::SimplicialLDLT<eigen_sparse, Eigen::Lower> eigen_ldlt;
   //
   const size_t    n_row_;         // number of rows (and columns) in H
   bool            init_done_;     // has init been called
   bool            update_called_; // has update been called
   sparse_rc       H_rc_;          // sparsity pattern for H
   //
   eigen_ldlt*     ptr_;           // eigens ldlt factorization
   //
public:
   // ----------------------------------------------------------------------
   // non-const functions
   //
   // constructor
   ldlt_eigen(size_t n_row);
   // destructor
   ~ldlt_eigen(void);
   // init
   void init(const sparse_rc& hes_rc);
   // update
   bool update(const CppAD::sparse_rcv<s_vector, v_vector>& hes_rcv);
   // ----------------------------------------------------------------------
   // const functions
   //
   // pattern
   const sparse_rc& pattern(void) const;
   //
   // split
   void split(
      eigen_sparse& L ,
      eigen_vector& D ,
      eigen_perm&   P
   ) const;
   //
   // logdet
   Double logdet(size_t& negative) const;
   // solve
   void solve_H(
      const CppAD::vector<size_t>& row     ,
      const CppAD::vector<Double>& val_in  ,
      CppAD::vector<Double>&       val_out
   ) const;
   //
   // sim_cov
   bool sim_cov(
   const CppAD::vector<Double>& w  ,
   CppAD::vector<Double>&       v
   ) const;
   //
   // compute a subset of the inverse
   void inv(
      const CppAD::vector<size_t>& row      ,
      const CppAD::vector<size_t>& col      ,
      CppAD::vector<Double>&       val
   ) const;
   // ----------------------------------------------------------------------
   // static functions
   //
   // solve_LDLT
   static eigen_vector solve_LDLT(
      const eigen_sparse& L  ,
      const eigen_vector& D  ,
      const eigen_perm&   P  ,
      const eigen_vector& b
   );
};


} } // END_CPPAD_MIXED_NAMESPACE


# endif
