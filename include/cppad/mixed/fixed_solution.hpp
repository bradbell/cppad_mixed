# ifndef CPPAD_MIXED_FIXED_SOLUTION_HPP
# define CPPAD_MIXED_FIXED_SOLUTION_HPP
// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin fixed_solution}
{xrst_spell
   lagrange
}

Optimal Solution Returned by optimize_fixed
###########################################

Syntax
******
``CppAD::mixed::fixed_solution`` *solution*

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE
   // END_PROTOTYPE
}

Convention
**********
If a Lagrange multiplier is non-zero (zero), the correspond constraint
is active (is not active) at the optimal solution.
The values specified below are as in the
:ref:`optimize_fixed@solution` return by ``optimize_fixed`` .

fixed_opt
*********
The size of this field is
:ref:`derived_ctor@n_fixed` .
It is the final value (optimal value found) for the fixed effects.

fixed_lag
*********
The size of this field is
:ref:`derived_ctor@n_fixed` .
If *solution* . ``fixed_lag`` [ *i* ]
is greater than zero (less than zero), it is
the Lagrange multiplier for the upper (lower) bound
for the *i*-th component of the fixed effects.

fix_con_lag
***********
The size of this field is
the number of fixed constraints; i.e., the size of
the vector :ref:`fix_constraint@vec` returned by
the ``fix_constraint`` function.
If *solution* . ``fix_con_lag`` [ *i* ]
is greater than zero (less than zero), it is
the Lagrange multiplier for the upper (lower) bound
for the *i*-th component of the fixed constraint function.

ran_con_lag
***********
The size of this field is
the number of random constraints; i.e.,
the number of rows
in the matrix :ref:`derived_ctor@A_rcv` .
If *solution* . ``ran_con_lag`` [ *i* ]
is greater than zero (less than zero), it is
the Lagrange multiplier for the upper (lower) bound
for the *i*-th row of the random constraint matrix :math:`A`.
{xrst_toc_hidden
   include/cppad/mixed/warm_start_struct.hpp
   include/cppad/mixed/trace_struct.hpp
}
warm_start
**********
This :ref:`warm_start_struct-name` contains
the necessary information to continue the ipopt optimization
from the current solution; i.e., warm start the optimization.

trace_vec
*********
The *i*-th element of this vector is a :ref:`trace_struct-name`
with the information corresponding to the *i*-th iteration
of the optimization algorithm.

{xrst_end fixed_solution}
------------------------------------------------------------------------------
*/
# include <cppad/utility/vector.hpp>
# include <cppad/mixed/warm_start_struct.hpp>
# include <cppad/mixed/trace_struct.hpp>

namespace CppAD { namespace mixed {
   // BEGIN_PROTOTYPE
   struct fixed_solution {
      CppAD::vector<double>       fixed_opt;
      CppAD::vector<double>       fixed_lag;
      CppAD::vector<double>       fix_con_lag;
      CppAD::vector<double>       ran_con_lag;
      warm_start_struct           warm_start;
      CppAD::vector<trace_struct> trace_vec;
   };
   // END_PROTOTYPE
} }


# endif
