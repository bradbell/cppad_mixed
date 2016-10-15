// $Id:$
# ifndef CPPAD_MIXED_FIXED_SOLUTION_HPP
# define CPPAD_MIXED_FIXED_SOLUTION_HPP
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
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
$$

$section Optimal Solution Returned by optimize_fixed$$

$head Syntax$$
$codei%CppAD::mixed::fixed_solution %solution%$$

$head Private$$
This structure is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Convention$$
If a Lagrange multiplier is zero, the correspond constraint
is not active at the optimal solution.

$head fixed_opt$$
This field has prototype
$codei%
	CppAD::vector<double> %solution%.fixed_opt
%$$
It has size zero when constructed
and when set its size is $cref/n_fixed/derived_ctor/n_fixed/$$.
It is the final value (optimal value found) for the fixed effects.

$head fixed_lag$$
This field has prototype
$codei%
	CppAD::vector<double> %solution%.fixed_lag
%$$
It has size zero when constructed
and when set its size is $cref/n_fixed/derived_ctor/n_fixed/$$.
If $icode%solution%.fixed_lag[%i%]%$$
is greater than zero (less than zero), it is
the Lagrange multiplier for the upper (lower) bound
for the $th i$$ component of the fixed effects.

$subhead Bounds$$
The size of the fixed effects upper and lower limits
$cref/fixed_upper/optimize_fixed/fixed_upper/$$ and
$cref/fixed_lower/optimize_fixed/fixed_lower/$$ are
used to determine the scale for the values of the fixed constraint function.
If for the $th i$$ component function,
one of the limits is infinite and the other is zero,
this scale cannot be determined and the values of
$icode%solution%.fix_[%i%]%$$ may not properly indicate
if the constraint is active.
In addition, you may get a warning about the solution check failing.

$head fix_con_lag$$
This field has prototype
$codei%
	CppAD::vector<double> %solution%.fix_con_lag
%$$
It has size zero when constructed
and when set its size is the number of fixed constraints; i.e., the size of
the vector $cref/vec/fix_constraint/vec/$$ returned by
the $code fix_constraint$$ function.
If $icode%solution%.fix_con_lag[%i%]%$$
is greater than zero (less than zero), it is
the Lagrange multiplier for the upper (lower) bound
for the $th i$$ component of the fixed constraint function.

$subhead Bounds$$
The size of the fixed constraints upper and lower limits
$cref/fix_constraint_upper/optimize_fixed/fix_constraint_upper/$$ and
$cref/fix_constraint_lower/optimize_fixed/fix_constraint_lower/$$ are
used to determine the scale for the values of the fixed constraint function.
If for the $th i$$ component function,
one of the limits is infinite and the other is zero,
this scale cannot be determined and the values of
$icode%solution%.fix_con_lag[%i%]%$$ may not properly indicate
if the constraint is active.
In addition, you may get a warning about the solution check failing.


$head ran_con_lag$$
This field has prototype
$codei%
	CppAD::vector<double> %solution%.ran_con_lag
%$$
It has size zero when constructed
and when set its size is the number of random constraints; i.e.,
the number of rows
in the matrix $cref/A_info/initialize/A_info/$$.
If $icode%solution%.ran_con_lag[%i%]%$$
is greater than zero (less than zero), it is
the Lagrange multiplier for the upper (lower) bound
for the $th i$$ row of the random constraint matrix $latex A$$.

$end
------------------------------------------------------------------------------
*/
# include <cppad/utility/vector.hpp>

namespace CppAD { namespace mixed {
	struct fixed_solution {
		CppAD::vector<double>  fixed_opt;
		CppAD::vector<double>  fixed_lag;
		CppAD::vector<double>  fix_con_lag;
		CppAD::vector<double>  ran_con_lag;
	};
} }


# endif
