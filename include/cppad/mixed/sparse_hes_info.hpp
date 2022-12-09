// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_SPARSE_HES_INFO_HPP
# define CPPAD_MIXED_SPARSE_HES_INFO_HPP
/*
{xrst_begin sparse_hes_info dev}
{xrst_spell
   work work
}

Sparse Hessian Information
##########################

Syntax
******
``CppAD::mixed::sparse_hes_info`` *hes_info*

Private
*******
This structure is an implementation detail and not part of the
CppAD Mixed user API.

Purpose
*******
This structure holds information about a specific Hessian.

row
***
The field *hes_info* . ``row`` has prototype

   ``CppAD::vector<size_t>`` *hes_info* . ``row``

It has size zero when it is constructed.
After initialization it should contain the row indices
for elements of the Hessian that get computed.

col
***
The field *hes_info* . ``col`` has prototype

   ``CppAD::vector<size_t>`` *hes_info* . ``col``

It has size zero when it is constructed.
After initialization it should contain the column indices
for elements of the Hessian that get computed.
It must has the same size as *hes_info* . ``row`` .

val
***
The field *hes_info* . ``val`` has prototype

   ``CppAD::vector<double>`` *hes_info* . ``val``

It has size zero when it is constructed.
After initialization it should be that same size as *hes_info* . ``row`` .
(In the special case where *hes_info* is used to compute
Hessian values with type other than ``double`` ,
*hes_info* . ``val`` can have size zero.)

work
****
The field *hes_info* . ``work`` has prototype

   ``CppAD::sparse_hessian_work work``

It has no information when it is constructed.
After initialization it should contain the CppAD cached information.

Sparse Hessian Call
*******************

| |tab| *f* . ``SparseHessian`` (
| |tab| |tab| *x* ,
| |tab| |tab| *w*
| |tab| |tab| *not_used* ,
| |tab| |tab| *hes_info* . ``row`` ,
| |tab| |tab| *hes_info* . ``col`` ,
| |tab| |tab| *hes_info* . ``val`` ,
| |tab| |tab| *hes_info* . ``work``
| |tab| )

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

not_used
========
This argument has one of the following prototypes

| |tab| ``CppAD::vector < std::set<size_t> >`` *not_used*
| |tab| ``CppAD::vector<bool>`` *not_used*
| |tab| ``CppAD::vectorBool`` *not_used*

It does not matter which and it is not used.

hes_info.val
============
Upon return *hes_info* . ``val`` [ *k* ] is the Hessian at row index
*hes_info* . ``row`` [ *k* ] and column index
*hes_info* . ``col`` [ *k* ] .

{xrst_end sparse_hes_info}
*/
namespace CppAD { namespace mixed {
   struct sparse_hes_info {
      CppAD::vector<size_t>      row;
      CppAD::vector<size_t>      col;
      CppAD::vector<double>      val;
      CppAD::sparse_hessian_work work;
   };
} }

# endif
