// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_SPARSE_JAC_RCV_HPP
# define CPPAD_MIXED_SPARSE_JAC_RCV_HPP
/*
{xrst_begin sparse_jac_rcv dev}
{xrst_spell
   nnz
}

Sparse Jacobian Computation Structure
#####################################

Syntax
******
``CppAD::mixed::sparse_jac_rcv`` *jac_rcv*

Private
*******
This structure is an implementation detail and not part of the
CppAD Mixed user API.

Purpose
*******
This structure holds information about a specific Jacobian

.. math::

   J(x) = f^{(1)} (x)

subset
******
The field *jac_rcv* . ``subset`` has prototype

   ``CppAD::mixed::d_sparse_rcv`` *jac_rcv* . ``subset``

It is empty zero when it is constructed; i.e.,
all of its sizes are zero.
After initialization it corresponds to the subset of the Jacobian
that is computed.

nnz
===
We use the notation

   *nnz* = *jac_rcv* . ``subset.nnz`` ()

row
===
We use the notation

   *row* = *jac_rcv* . ``subset.row`` ()

col
===
We use the notation

   *col* = *jac_rcv* . ``subset.col`` ()

val
===
We use the notation

   *val* = *jac_rcv* . ``subset.val`` ()

forward
*******
The field *jac_rcv* . ``forward`` has prototype

   *bool* ``jac_rcv`` . *forward*

work
****
The field *jac_rcv* . ``work`` has prototype

   ``CppAD::sparse_jac_work`` *jac_rcv* . ``work``

It has no information when it is constructed; i.e., it is empty.
After initialization it should contain the CppAD cached information
for the call to ``sparse_jac_for`` or ``sparse_jac_rev``
specified below.

Computing Sparse Jacobians
**************************
Upon return (from ``sparse_jac_for`` or ``sparse_jac_rev`` ),
for *k* = 0 , ..., *nnz* ``-1`` ,
*val* [ *k* ] is the value of :math:`J_{i,j} (x)`
where *i* = *row* [ *k* ] and *j* = *col* [ ``k`` ] .

Forward Mode
============
If *jac_rcv* . ``forward`` is true,

| |tab| *f* . ``sparse_jac_for`` (
| |tab| |tab| *group_max* ,
| |tab| |tab| *x* ,
| |tab| |tab| *jac_rcv* . ``subset`` ,
| |tab| |tab| *not_used_pattern* ,
| |tab| |tab| *not_used_coloring* ,
| |tab| |tab| *jac_rcv* . ``work``
| |tab| )

computes the Jacobian values.

Reverse Mode
============
If *jac_rcv* . ``forward`` is false,

| |tab| *f* . ``sparse_jac_for`` (
| |tab| |tab| *x* ,
| |tab| |tab| *jac_rcv* . ``subset`` ,
| |tab| |tab| *not_used_pattern* ,
| |tab| |tab| *not_used_coloring* ,
| |tab| |tab| *jac_rcv* . ``work``
| |tab| )

computes the Jacobian values.

f
=
The function *f* has prototype

   ``CppAD::ADFun<double>`` *f*

group_max
=========
The parameter *group_max* has prototype

   ``size_t`` *group_max*

and must be greater than zero.
It specifies the maximum number of directions to group during
a single forward sweep.
This uses separate memory for each direction (more memory),
but my be significantly faster.

x
=
The argument *x* has prototype

   ``const CppAD::vector<double>&`` *x*

and its size is *f* . ``Domain`` () .
It is the location where the Jacobian is being evaluated.

not_used_pattern
================
This argument has the following prototype

   ``const CppAD::mixed::sparse_rc&`` *not_used_pattern*

It is not used and hence its value does not matter.

not_used_coloring
=================
This argument has the following prototype

   ``const std::string&`` *not_used_coloring*

It is not used and hence its value does not matter.

{xrst_end sparse_jac_rcv}
*/
# include <cppad/mixed/typedef.hpp>

namespace CppAD { namespace mixed {
   struct sparse_jac_rcv {
      bool                   forward;
      d_sparse_rcv           subset;
      CppAD::sparse_jac_work work;
   };
} }

# endif
