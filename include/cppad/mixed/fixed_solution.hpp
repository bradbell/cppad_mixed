// $Id:$
# ifndef CPPAD_MIXED_FIXED_SOLUTION_HPP
# define CPPAD_MIXED_FIXED_SOLUTION_HPP
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */

/*
$begin fixed_solution$$
$spell
	CppAD
	vec
	rcv
	ipopt
$$

$section Optimal Solution Returned by optimize_fixed$$

$head Syntax$$
$codei%CppAD::mixed::fixed_solution %solution%$$

$head Prototype$$
$srcthisfile%0%// BEGIN_PROTOTYPE%// END_PROTOTYPE%1%$$

$head Public$$
This structure is part of the
CppAD Mixed user API.

$head Convention$$
If a Lagrange multiplier is non-zero (zero), the correspond constraint
is active (is not active) at the optimal solution.
The values specified below are as in the
$cref/solution/optimize_fixed/solution/$$ return by $code optimize_fixed$$.

$head fixed_opt$$
The size of this field is
$cref/n_fixed/derived_ctor/n_fixed/$$.
It is the final value (optimal value found) for the fixed effects.

$head fixed_lag$$
The size of this field is
$cref/n_fixed/derived_ctor/n_fixed/$$.
If $icode%solution%.fixed_lag[%i%]%$$
is greater than zero (less than zero), it is
the Lagrange multiplier for the upper (lower) bound
for the $th i$$ component of the fixed effects.

$head fix_con_lag$$
The size of this field is
the number of fixed constraints; i.e., the size of
the vector $cref/vec/fix_constraint/vec/$$ returned by
the $code fix_constraint$$ function.
If $icode%solution%.fix_con_lag[%i%]%$$
is greater than zero (less than zero), it is
the Lagrange multiplier for the upper (lower) bound
for the $th i$$ component of the fixed constraint function.

$head ran_con_lag$$
The size of this field is
the number of random constraints; i.e.,
the number of rows
in the matrix $cref/A_rcv/derived_ctor/A_rcv/$$.
If $icode%solution%.ran_con_lag[%i%]%$$
is greater than zero (less than zero), it is
the Lagrange multiplier for the upper (lower) bound
for the $th i$$ row of the random constraint matrix $latex A$$.

$children%
	include/cppad/mixed/warm_start_struct.hpp
%$$
$head warm_start$$
This $cref warm_start_struct$$ contains
the necessary information to continue the ipopt optimization
from the current solution; i.e., warm start the optimization.

$end
------------------------------------------------------------------------------
*/
# include <cppad/utility/vector.hpp>
# include <cppad/mixed/warm_start_struct.hpp>

namespace CppAD { namespace mixed {
	// BEGIN_PROTOTYPE
	struct fixed_solution {
		CppAD::vector<double>  fixed_opt;
		CppAD::vector<double>  fixed_lag;
		CppAD::vector<double>  fix_con_lag;
		CppAD::vector<double>  ran_con_lag;
		warm_start_struct      warm_start;
	};
	// END_PROTOTYPE
} }


# endif
