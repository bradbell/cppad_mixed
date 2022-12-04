/*
// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
{xrst_begin information_mat}

Compute the Observed Information For Fixed Effects
##################################################

Deprecated 2020-03-22
*********************
Use :ref:`hes_fixed_obj-name` instead.

Syntax
******

| *information_rcv* = *mixed_object* . ``information_mat`` (
| |tab| *solution* , *random_opt*
| )

Purpose
*******
Compute the observed information matrix.
We use :math:`L ( \theta )` to denote the
:ref:`fixed effects objective<theory@Objective@Fixed Effects Objective, L(theta)>` .
The observed information is

.. math::

   L^{(2)} ( \hat{\theta} )

Absolute value terms in the
:ref:`problem@Negative Log-Density Vector`
for the :ref:`fix_likelihood-name` are not include in this Hessian
(because they do not have a derivative, let alone Hessian, at zero).

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

solution
********
is the :ref:`optimize_fixed@solution`
for a previous call to :ref:`optimize_fixed-name` .
Only the *solution* . ``fixed_opt`` field is used.

random_opt
**********
is the optimal random effects corresponding to the solution; i.e.

| |tab| *random_opt* = *mixed_object* . ``optimize_random`` (
| |tab| |tab| *random_options* ,
| |tab| |tab| *solution* . ``fixed_opt`` ,
| |tab| |tab| *random_lower* ,
| |tab| |tab| *random_upper* ,
| |tab| |tab| *random_in*
| |tab| )

*random_options* ,
*random_lower* ,
*random_upper* , and
*random_in* , are the same
as in the call to ``optimize_fixed`` that corresponds to *solution* .

information_rcv
***************
The return value has prototype

   ``CppAD::mixed::d_sparse_rcv`` *information_rcv*

see :ref:`typedef@Sparse Types@d_sparse_rcv` .
This is a sparse matrix representation for the
lower triangle of the observed information matrix,
which is symmetric and hence determined by its lower triangle.
Absolute value terms in the
:ref:`problem@Negative Log-Density Vector`
for the :ref:`fix_likelihood-name` are not include in this Hessian
because they do not have a derivative (let alone Hessian) at zero.
{xrst_toc_hidden
   example/user/information_mat.cpp
}

Example
*******
The file :ref:`information_mat.cpp-name` contains an example and
test of this routine. It returns true for success and false for failure.

{xrst_end information_mat}
------------------------------------------------------------------------------
*/
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/triple2eigen.hpp>
# include <cppad/mixed/exception.hpp>

// ---------------------------------------------------------------------------
CppAD::mixed::d_sparse_rcv cppad_mixed::information_mat(
   const CppAD::mixed::fixed_solution&  solution             ,
   const d_vector&                      random_opt           )
{  d_sparse_rcv result;
   try
   {  result = try_hes_fixed_obj(solution.fixed_opt, random_opt);
   }
   catch(const std::exception& e)
   {  std::string error_message = "information_mat: std::exception: ";
      error_message += e.what();
      fatal_error(error_message);
      assert(false);
   }
   catch(const CppAD::mixed::exception& e)
   {  std::string error_message = e.message("information_mat");
      fatal_error(error_message);
      assert(false);
   }
   //
   return result;
}
