// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_SPARSE_HES_RCV_HPP
# define CPPAD_MIXED_SPARSE_HES_RCV_HPP
/*
{xrst_begin sparse_hes_rcv dev}
{xrst_spell
  nnz
}

Sparse Hessian Computation Structure
####################################

Syntax
******
``CppAD::mixed::sparse_hes_rcv`` *hes_rcv*

Private
*******
This structure is an implementation detail and not part of the
CppAD Mixed user API.

Purpose
*******
This structure holds information about the Hessian

.. math::

   H(x) = \sum_{i=0}^{m-1} w_i f_i^{(2)} (x)

where :math:`m` is the number of components in the function :math:`f(x)`.

subset
******
The field *hes_rcv* . ``subset`` has prototype

   ``CppAD::mixed::d_sparse_rcv`` *hes_rcv* . ``subset``

It is empty zero when it is constructed; i.e.,
all of its sizes are zero.
After initialization it corresponds to the subset of the Hessian
that is computed.

nnz
===
We use the notation

   *nnz* = *hes_rcv* . ``subset.nnz`` ()

row
===
We use the notation

   *row* = *hes_rcv* . ``subset.row`` ()

col
===
We use the notation

   *col* = *hes_rcv* . ``subset.col`` ()

val
===
We use the notation

   *val* = *hes_rcv* . ``subset.val`` ()

work
****
The field *hes_rcv* . ``work`` has prototype

   ``CppAD::sparse_hes_work`` *hes_rcv* . ``work``

It has no information when it is constructed; i.e., it is empty.
After initialization it should contain the CppAD cached information
for the call to ``sparse_hes`` specified below.

Computing Sparse Hessians
*************************

| |tab| *f* . ``sparse_hes`` (
| |tab| |tab| *x* ,
| |tab| |tab| *w*
| |tab| |tab| *hes_rcv* . ``subset`` ,
| |tab| |tab| *not_used_pattern* ,
| |tab| |tab| *not_used_coloring* ,
| |tab| |tab| *hes_rcv* . ``work``
| |tab| )

Upon return,
for *k* = 0 , ..., *nnz* ``-1`` ,
*val* [ *k* ] is the value of :math:`H_{i,j} (x)`
where *i* = *row* [ *k* ] and *j* = *col* [ ``k`` ] .

f
=
The function *f* has prototype

   ``CppAD::ADFun<double>`` *f*

x
=
The argument *x* has prototype

   ``const CppAD::vector<double>&`` *x*

and its size is *f* . ``Domain`` () .
It is the location where the Hessian is being evaluated.

w
=
The argument *w* has prototype

   ``const CppAD::vector<double>&`` *w*

and its size is *f* . ``Range`` () .
It is the weighting for the components of the Hessian that is
being computed.

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

{xrst_end sparse_hes_rcv}
*/
# include <cppad/mixed/typedef.hpp>

namespace CppAD { namespace mixed {
   struct sparse_hes_rcv {
      d_sparse_rcv           subset;
      CppAD::sparse_hes_work work;
   };
} }

# endif
