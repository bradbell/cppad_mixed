# ifndef CPPAD_MIXED_WARM_START_STRUCT_HPP
# define CPPAD_MIXED_WARM_START_STRUCT_HPP
// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------

/*
$begin warm_start_struct$$
$spell
   mu
   Ipopt
   struct
   cppad
   CppAD
$$

$section Ipopt Warm Start Information$$

$head Syntax$$
$codei%CppAD::mixed::warm_start_struct %warm_start%$$

$head Prototype$$
$srcthisfile%0%// BEGIN_PROTOTYPE%// END_PROTOTYPE%1%$$

$head mu$$
This is the warm start value for the barrier penalty parameter.

$head scale_f$$
This is the cppad_mixed scaling factor for the ipopt objective function
$latex f(x)$$.

$head x_info$$
If the size of this vector is zero, the size of $icode g_info$$
must also be zero.
Otherwise, $code x_info$$ has size equal to the number of primal variables
$latex x$$.
The $th j$$ element of this vector contains the following fields:

$subhead x$$
is the warm start value for x[j].

$subhead z_L$$
is the warm start value for z_L[j].

$subhead z_U$$
is the warm start value for z_U[j].

$head scale_x$$
is the cppad_mixed scaling factor for x[j].

$head g_info$$
If the size of $icode x_info$$ is non-zero,
$code g_info$$ has size equal to the number of $latex g(x)$$ constraints.
The $th i$$ element of this vector contains the following fields:

$subhead lambda$$
is the warm start value for lambda[i].

$subhead scale_g$$
is the cppad_mixed scaling factor for $latex g_i(x)$$.

$head Public$$
This structure is part of the CppAD Mixed user API.

$end
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
