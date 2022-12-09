// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/ipopt_fixed.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
{xrst_begin ipopt_fixed_new_random dev}

Compute New Random Effects and Update Factor
############################################

Syntax
******
``new_random`` ( *fixed_vec* )

ipopt_fixed
***********
This is a private member function of the :ref:`ipopt_fixed-name` class.

n_random\_
**********
Is assumed that this member variable is greater than zero.

random_ipopt_options\_
**********************
This member variable contains
the value of the
:ref:`ipopt_fixed_ctor@random_ipopt_options`
in the ``ipopt_fixed`` constructor.

random_lower\_
**************
This member variable contains
the value of the :ref:`ipopt_fixed_ctor@random_lower`
in the ``ipopt_fixed`` constructor.

random_upper\_
**************
This member variable contains
the value of the :ref:`ipopt_fixed_ctor@random_upper`
in the ``ipopt_fixed`` constructor.

random_in\_
***********
This member variable contains
the value of the :ref:`ipopt_fixed_ctor@random_in`
in the ``ipopt_fixed`` constructor.

fixed_vec
*********
This argument has prototype

   ``const d_vector&`` *fixed_vec*

it is the value of the fixed effects that we are computing the random
effects and updated factor for.

random_cur\_
************
This member variable contain is set the optimal random effects
corresponding to *fixed_vec* .

mixed_object\_
**************
The factor in this member variables is updated using the call

   ``mixed_object_.update_factor`` ( *fixed_vec* , ``random_cur_`` )

see :ref:`update_factor-name` for side effects.

Prototype
*********
{xrst_spell_off}
{xrst_code cpp} */
void ipopt_fixed::new_random(const d_vector& fixed_vec)
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_fixed_new_random}
*/
{  assert( n_random_ > 0 );
   // Compute the optimal random effects corresponding to fixed effects.
   // Use try_optimize_random instead of optimize_random, so that thows
   // are caught at the fixed effects, not random effects, level.
   random_cur_ = mixed_object_.try_optimize_random(
      random_ipopt_options_,
      fixed_vec,
      random_lower_,
      random_upper_,
      random_in_
   );
   mixed_object_.update_factor(fixed_vec, random_cur_);
}
} } // END_CPPAD_MIXED_NAMESPACE
