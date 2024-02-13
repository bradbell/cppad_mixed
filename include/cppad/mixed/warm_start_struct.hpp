# ifndef CPPAD_MIXED_WARM_START_STRUCT_HPP
# define CPPAD_MIXED_WARM_START_STRUCT_HPP
// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------

/*
{xrst_begin warm_start_struct}

Ipopt Warm Start Information
############################

Syntax
******
``CppAD::mixed::warm_start_struct`` *warm_start*

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE
   // END_PROTOTYPE
}

mu
**
This is the warm start value for the barrier penalty parameter.

scale_f
*******
This is the cppad_mixed scaling factor for the ipopt objective function
:math:`f(x)`.

x_info
******
If the size of this vector is zero, the size of *g_info*
must also be zero.
Otherwise, ``x_info`` has size equal to the number of primal variables
:math:`x`.
The *j*-th element of this vector contains the following fields:

x
=
is the warm start value for x[j].

z_L
===
is the warm start value for z_L[j].

z_U
===
is the warm start value for z_U[j].

scale_x
*******
is the cppad_mixed scaling factor for x[j].

g_info
******
If the size of *x_info* is non-zero,
``g_info`` has size equal to the number of :math:`g(x)` constraints.
The *i*-th element of this vector contains the following fields:

lambda
======
is the warm start value for lambda[i].

scale_g
=======
is the cppad_mixed scaling factor for :math:`g_i(x)`.

Public
******
This structure is part of the CppAD Mixed user API.

{xrst_end warm_start_struct}
------------------------------------------------------------------------------
*/
# include <cppad/utility/vector.hpp>

namespace CppAD { namespace mixed {
   // BEGIN_PROTOTYPE
   struct x_info_struct {
      double x; double z_L; double z_U; double scale_x;
   };
   struct g_info_struct {
      double lambda; double scale_g;
   };
   struct warm_start_struct {
      double mu;
      double scale_f;
      CppAD::vector<x_info_struct> x_info;
      CppAD::vector<g_info_struct> g_info;
   };
   // END_PROTOTYPE
} }


# endif
