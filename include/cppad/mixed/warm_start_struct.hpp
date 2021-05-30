// $Id:$
# ifndef CPPAD_MIXED_WARM_START_STRUCT_HPP
# define CPPAD_MIXED_WARM_START_STRUCT_HPP
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */

/*
$begin warm_start_struct$$
$spell
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

$head scale_f$$
This is the cppad_mixed scaling factor for the ipopt objective function
$latex f(x)$$.

$head x_info$$
This vector has size equal to the number of primal variables $latex x$$.
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
This vector has size equal to the number of $latex g(x)$$ constraints.
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
		double scale_f;
		CppAD::vector<x_info_struct> x_info;
		CppAD::vector<g_info_struct> g_info;
	};
	// END_PROTOTYPE
} }


# endif
