// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin ldlt_cholmod_pattern dev}

Update Factorization Using new Matrix Values
############################################

Syntax
******
*H_rc* = *ldlt_obj* . ``pattern`` ()

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE
   // END_PROTOTYPE
}

Private
*******
The :ref:`ldlt_cholmod-name` class is an
:ref:`implementation detail<ldlt_cholmod@Private>` and not part of the
CppAD Mixed user API.

ldlt_obj
********
This object has prototype

   ``CppAD::mixed::ldlt_cholmod`` *ldlt_obj*

In addition, it must have a previous call to
:ref:`ldlt_cholmod_init-name` .

H_rc
****
The return value is a copy of the sparsity pattern
:ref:`ldlt_cholmod_init@H_rc` in the corresponding call to
*ldlt_obj* . ``init`` ( *H_rc* ) .

Order of Operations
*******************
This *ldlt_obj* function must be called,
after the constructor and :ref:`init<ldlt_cholmod_init-name>` .

Example
*******
The file :ref:`ldlt_cholmod.cpp<ldlt_cholmod.cpp@pattern>` contains an
example and test that uses this function.

{xrst_end ldlt_cholmod_pattern}
*/
// ----------------------------------------------------------------------------


# include <cppad/mixed/include_cholmod.hpp>
# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cassert>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN_PROTOTYPE
const sparse_rc& ldlt_cholmod::pattern(void) const
// END_PROTOTYPE
{  assert( init_done_ );
   return H_rc_;
}

} } // END_CPPAD_MIXED_NAMESPACE
