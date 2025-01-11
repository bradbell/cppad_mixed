// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin ran_con_jac dev}
{xrst_spell
  nnz
}

Jacobian of the Random Constraint Function
##########################################

Syntax
******

| *mixed_object* . ``ran_con_jac`` (
| |tab| *fixed_vec* , *random_vec* , *jac_rcv*
| )

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

fixed_vec
*********
This argument has prototype

   ``const CppAD::vector<double>&`` *fixed_vec*

It specifies the value of the
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>`
vector :math:`\theta`.

random_vec
**********
This argument has prototype

   ``const CppAD::vector<double>&`` *random_vec*

It specifies the value of the
:ref:`random effects<problem@Notation@Random Effects, u>`
vector :math:`u`,
This must be the
:ref:`optimal random effects<theory@Optimal Random Effects, u^(theta)>`
:math:`\hat{u} ( \theta )`.

jac_rcv
*******
This argument has prototype

   ``d_sparse_rcv&`` *jac_rcv*

see :ref:`typedef@Sparse Types@d_sparse_rcv` .

Sparsity Pattern
================
In the case where the input value of *jac_rcv* is empty
(all its sizes are zero)
only the :ref:`typedef@Sparse Types@sparse_rc`
sparsity pattern in *jac_rcv* computes.
To be specific, upon return *jac_rcv* . ``val`` () has size
equal to the number of non-zero elements in the sparsity pattern,
but the value of the elements is not specified.

Sparse Matrix
=============
In the case where the input value of *jac_rcv* is non-empty
upon return *jac_rcv* is a sparse matrix representation
of the Jacobian of the
:ref:`random constraint function<problem@Notation@Random Constraint Function, A*u^(theta)>` ; i.e.,

.. math::

   \partial_\theta [ A \hat{u} ( \theta ) ]
   =
   A \hat{u}_\theta ( \theta )
   =
   - A
   f_{u,u} \left[ \theta , \hat{u} ( \theta ) \right]^{-1}
   f_{u,\theta} \left[ \theta , \hat{u} ( \theta )  \right]

see
:ref:`theory@Derivative of Optimal Random Effects` .

Column Major Order
==================
The results in *jac_rcv* are in column major order; i.e.,
upon return

   *jac_rcv* . ``col_major`` [ *k* ] = *k*

for *k* = 0, ..., *jac_rcv* . ``nnz`` () ``-1`` .
{xrst_toc_hidden
   example/private/ran_con_jac.cpp
}

Example
*******
The file :ref:`ran_con_jac.cpp-name` contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

{xrst_end ran_con_jac}
*/
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>

void cppad_mixed::ran_con_jac(
   const d_vector&                fixed_vec  ,
   const d_vector&                random_vec ,
   d_sparse_rcv&                  jac_rcv    )
{  assert( fixed_vec.size()  == n_fixed_ );
   assert( random_vec.size() == n_random_ );
   assert( A_rcv_.nr() != 0 ); // could just return here in this case
   //
   bool empty =  jac_rcv.nr()==0 && jac_rcv.nc()==0 && jac_rcv.nnz()==0;
   // 2DO: convert jac_rcv from a dense matrix to a truely sparse matrix
   if( empty )
   {  size_t nr  = A_rcv_.nr();
      size_t nc  = n_fixed_;
      size_t nnz = nr * nc;
      sparse_rc pattern(nr, nc, nnz);
      for(size_t j = 0; j < nc; j++)
      {  for(size_t i = 0; i < nr; i++)
            pattern.set(j * nr + i, i, j);
      }
      //
      jac_rcv = d_sparse_rcv(pattern);
      //
      return;
   }
   //
   // packed version of fixed and random effects
   d_vector both(n_fixed_ + n_random_);
   pack(fixed_vec, random_vec, both);

   // Compute the Hessian cross terms f_{u,theta} ( theta , u )
   sparse_rc  not_used_pattern;
   std::string not_used_coloring;
   size_t K = hes_cross_.subset.nnz();
   d_vector w(1);
   w[0] = 1.0;
   ran_like_fun_.sparse_hes(
      both,
      w,
      hes_cross_.subset,
      not_used_pattern,
      not_used_coloring,
      hes_cross_.work
   );
   if( CppAD::hasnan( hes_cross_.subset.val() ) ) CppAD::mixed::exception(
      "ran_con_jac", "result has a nan"
   );
   //
   // Use the column major order specification for
   // (hes_cross_.row, hes_cross_.col)
   size_t k = 0;
   size_t row = n_random_;
   size_t col = n_fixed_;
   if( k < K )
   {  assert( hes_cross_.subset.row()[k] >= n_fixed_ );
      row = hes_cross_.subset.row()[k] - n_fixed_;
      col = hes_cross_.subset.col()[k];
      assert( row < n_random_ );
      assert( col < n_fixed_ );
   }
   //
   // 2DO: figure out how to make this a sparse calculation.
   // The problem is propagating the sparsity though Choleksy inverse.
   CppAD::vector<size_t> row_solve(n_random_);
   d_vector val_b(n_random_), val_x(n_random_), Au_theta_j(A_rcv_.nr());
   for(size_t i = 0; i < n_random_; i++)
   {  row_solve[i] = i;
   }
   // Loop over fixed effects
   size_t ell = 0;
   for(size_t j = 0; j < n_fixed_; j++)
   {  // b = j-th column of - f_{u, theta} (theta, u)
      for(size_t i = 0; i < n_random_; i++)
         val_b[i]     = 0.0;
      while( col <= j )
      {  assert( col == j );
         val_b[row] = -hes_cross_.subset.val()[k];
         //
         k++;
         if( k < K )
         {  assert( hes_cross_.subset.row()[k] >= n_fixed_ );
            row = hes_cross_.subset.row()[k] - n_fixed_;
            col = hes_cross_.subset.col()[k];
            assert( row < n_random_ );
            assert( col < n_fixed_ );
         }
         else
         {  row = n_random_;
            col = n_fixed_;
         }
      }
      assert( col > j );
      // x = j-th column of - f_{u,u}(theta, u)^{-1} f_{u,theta}(theta, u)
      //   = j-th column of uhat_theta ( theta )
      ldlt_ran_hes_.solve_H(row_solve, val_b, val_x);
      //
      // multipliy A times j-th column of uhat_theta
      for(size_t i = 0; i < A_rcv_.nr(); i++)
         Au_theta_j[i] = 0.0;
      size_t M = A_rcv_.row().size();
      for(size_t m = 0; m < M; m++)
      {  size_t r = A_rcv_.row()[m];
         size_t c = A_rcv_.col()[m];
         double v = A_rcv_.val()[m];
         Au_theta_j[r] += v * val_x[c];
      }
      for(size_t i = 0; i < A_rcv_.nr(); i++)
      {  // computed values
         assert( jac_rcv.row()[ell] == i );
         assert( jac_rcv.col()[ell] == j );
         jac_rcv.set(ell, Au_theta_j[i]);
         ++ell;
      }
   }
   return;
}
