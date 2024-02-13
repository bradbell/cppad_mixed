# ifndef CPPAD_MIXED_TRACE_STRUCT_HPP
# define CPPAD_MIXED_TRACE_STRUCT_HPP
// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------

/*
{xrst_begin trace_struct}
{xrst_spell
   du
   iter
   pr
}

Ipopt Trace Information
#######################

Syntax
******
``CppAD::mixed::trace_struct`` *trace*

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE
   // END_PROTOTYPE
}

iter
****
See ipopt trace :ref:`ipopt_trace@iter` .

obj_value
*********
See ipopt trace :ref:`ipopt_trace@objective` .

inf_pr
******
See ipopt trace :ref:`ipopt_trace@inf_pr` .

inf_du
******
See ipopt trace :ref:`ipopt_trace@inf_du` .

mu
**
See ipopt trace :ref:`ipopt_trace@lg(mu)` .

d_norm
******
See ipopt trace :ref:`ipopt_trace@||d||` .

regularization_size
*******************
See ipopt trace :ref:`ipopt_trace@lg(rg)` .

alpha_du
********
See ipopt trace :ref:`ipopt_trace@alpha_du` .

alpha_pr
********
See ipopt trace :ref:`ipopt_trace@alpha_pr` .

ls_trials
*********
See ipopt trace :ref:`ipopt_trace@ls` .

restoration
***********
Is ipopt currently in restoration mode.

{xrst_end trace_struct}
------------------------------------------------------------------------------
*/
namespace CppAD { namespace mixed {
   // BEGIN_PROTOTYPE
   struct trace_struct {
      size_t iter;
      double obj_value;
      double inf_pr;
      double inf_du;
      double mu;
      double d_norm;
      double regularization_size;
      double alpha_du;
      double alpha_pr;
      size_t ls_trials;
      bool   restoration;
   };
   // END_PROTOTYPE
} }


# endif
