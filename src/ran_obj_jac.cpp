// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>
/*
{xrst_begin ran_obj_jac}

Derivative of Laplace Objective
###############################

Syntax
******
*mixed_object* . ``ran_obj_jac`` ( *fixed_vec* , *random_vec* , *r_fixed* )

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

Purpose
*******
This routine computes the
:ref:`derivative of the Laplace objective<theory@Derivative of Laplace Objective>` .

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

ldlt_ran_hes\_
**************
It is assumed that the member variable

   ``CPPAD_MIXED_LDLT_CLASS ldlt_ran_hes_``

was updated using :ref:`update_factor-name` for the specified values of the
fixed and random effects.

fixed_vec
*********
This argument has prototype

   ``const CppAD::vector<double>&`` *fixed_vec*

and is the value of fixed effects :math:`\theta`.

random_vec
**********
This argument has prototype

   ``const CppAD::vector<double>&`` *random_vec*

and is the value of fixed effects :math:`u`.
It is assumed that these random effects are optimal for the specified
fixed effects and hence :math:`f_u ( \theta , u ) = 0`.

r_fixed
*******
This argument has prototype

   ``CppAD::vector<double>&`` *r_fixed*

If the input size must be equal to ``n_fixed_`` .
Upon return, it contains the value of the derivative w.r.t
the fixed effects; i.e. :math:`r_\theta ( \theta )`.
{xrst_toc_hidden
   example/private/ran_obj_jac.cpp
}
Example
*******
The file :ref:`ran_obj_jac.cpp-name` contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

{xrst_end ran_obj_jac}
*/
// ----------------------------------------------------------------------------
void cppad_mixed::ran_obj_jac(
   const d_vector& fixed_vec  ,
   const d_vector& random_vec ,
   d_vector&       r_fixed    )
{
   assert( fixed_vec.size() == n_fixed_ );
   assert( random_vec.size() == n_random_ );
   assert( r_fixed.size() == n_fixed_ );

   //
   // Compute derivative of logdet of f_{u,u} ( theta , u )
   d_vector logdet_fix(n_fixed_), logdet_ran(n_random_);
   logdet_jac(fixed_vec, random_vec, logdet_fix, logdet_ran);

   // packed version of fixed and random effects
   d_vector both(n_fixed_ + n_random_);
   pack(fixed_vec, random_vec, both);

   //
   // Compute derivative of f(theta , u ) w.r.t theta and u
   d_vector w(1), f_both(n_fixed_ + n_random_);
   w[0] = 1.0;
   ran_like_fun_.Forward(0, both);
   f_both = ran_like_fun_.Reverse(1, w);
   if( CppAD::hasnan( f_both ) ) CppAD::mixed::exception(
      "ran_obj_jac", "result has a nan"
   );
   d_vector f_fixed(n_fixed_), f_random(n_random_);
   unpack(f_fixed, f_random, f_both);
   //
   // Compute the Hessian cross terms f_{u,theta} ( theta , u )
   // 2DO: another ran_like_fun_.Forward(0, both) is done by SparseHessian
   std::string not_used_coloring;
   sparse_rc   not_used_pattern;
   size_t K = hes_cross_.subset.nnz();
   ran_like_fun_.sparse_hes(
      both,
      w,
      hes_cross_.subset,
      not_used_pattern,
      not_used_coloring,
      hes_cross_.work
   );
   if( CppAD::hasnan( hes_cross_.subset.val() ) ) CppAD::mixed::exception(
      "ran_obj_jac", "result has a nan"
   );
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

   // Loop over fixed effects and compute r_fixed one component at a time
   CppAD::vector<bool> found(n_random_);
   CppAD::vector<size_t> row_solve;
   CppAD::vector<double> val_b, val_x;
   for(size_t j = 0; j < n_fixed_; j++)
   {  // vectors for this column
      row_solve.resize(0);
      val_b.resize(0);
      val_x.resize(0);
      //
      // increment through all row indices
      size_t i = 0;

      // b = j-th column of - f_{u, theta} (theta, u)
      while( col <= j )
      {  // make sure we are solving for all rows where
         // derivative of logdet w.r.t random effects is non-zero.
         while( i < row )
         {  if( logdet_ran[i] != 0.0 )
            {  row_solve.push_back(i);
               val_b.push_back(0.0);
            }
            i++;
         }
         if( col == j )
         {  row_solve.push_back(row);
            val_b.push_back( - hes_cross_.subset.val()[k] );
            found[row] = true;
            assert( i == row );
            i++;
         }
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
      // make sure we are solving for all rows where
      // derivative of logdet w.r.t random effects is non-zero.
      while( i < n_random_ )
      {  if( logdet_ran[i] != 0.0 )
         {  row_solve.push_back(i);
            val_b.push_back(0.0);
         }
         i++;
      }

      // x = j-th column of - f_{u,u}(theta, u)^{-1} f_{u,theta}(theta, u)
      val_x.resize( row_solve.size() );
      ldlt_ran_hes_.solve_H(row_solve, val_b, val_x);
      //
      // parial w.r.t fixed effects contribution to total derivative
      r_fixed[j] =  f_fixed[j] + 0.5 * logdet_fix[j];
      //
      // compute effect of uhat_{theta(j)} (theta) on the total derivative
      for(size_t ell = 0; ell < row_solve.size(); ell++)
      {  // random effect index
         i = row_solve[ell];
         // partial of optimal random effect for this (i, j)
         double ui_thetaj = val_x[ell];
         // partial of random part of objective w.r.t this random effect
         // Note f_random = 0 because u is optimal for the fixed effects
         // double h_ui = f_random[i] + 0.5 * logdet_ran[i];
         double h_ui = 0.5 * logdet_ran[i];
         // contribution to total derivative
         r_fixed[j] += h_ui * ui_thetaj;
      }
   }
   //
   return;
}
