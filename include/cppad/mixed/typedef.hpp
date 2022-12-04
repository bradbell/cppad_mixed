// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_TYPEDEF_HPP
# define CPPAD_MIXED_TYPEDEF_HPP

# include <cppad/cppad.hpp>

/*
{xrst_begin typedef}
{xrst_spell
   hpp
   namespace
   typedef
}

Types Defined in the CppAD Mixed Namespace
##########################################

Syntax
******

   # ``include <cppad/mixed/typedef.hpp>``

Begin Namespace
***************
All the definitions below are made inside the ``CppAD::mixed`` namespace;
i.e.,
{xrst_spell_off}
{xrst_code cpp} */
namespace CppAD { namespace mixed {

/* {xrst_code}
{xrst_spell_on}

Scalar Types
************

a1_double
=========
Scalar with one level of AD:
{xrst_spell_off}
{xrst_code cpp} */
   typedef CppAD::AD<double> a1_double;
/* {xrst_code}
{xrst_spell_on}

Vector Types
************

s_vector
========
Vectors with elements of type ``size_t`` :
{xrst_spell_off}
{xrst_code cpp} */
   typedef CppAD::vector<size_t> s_vector;
/* {xrst_code}
{xrst_spell_on}

d_vector
========
Vectors with elements of type ``double`` :
{xrst_spell_off}
{xrst_code cpp} */
   typedef CppAD::vector<double> d_vector;
/* {xrst_code}
{xrst_spell_on}

a1_vector
=========
Vectors with elements of that have one level of AD:
{xrst_spell_off}
{xrst_code cpp} */
   typedef CppAD::vector<a1_double> a1_vector;
/* {xrst_code}
{xrst_spell_on}

Sparse Types
************

sparse_rc
=========
Sparsity patterns using index vector of type ``s_vector`` :
{xrst_spell_off}
{xrst_code cpp} */
   typedef CppAD::sparse_rc<s_vector> sparse_rc;
/* {xrst_code}
{xrst_spell_on}

d_sparse_rcv
============
Sparse matrices using index vector of type ``s_vector``
and value vectors of type ``d_vector`` :
{xrst_spell_off}
{xrst_code cpp} */
   typedef CppAD::sparse_rcv<s_vector, d_vector> d_sparse_rcv;

/* {xrst_code}
{xrst_spell_on}
a1_sparse_rcv
=============
Sparse matrices using index vector of type ``s_vector``
and value vectors of type ``a1_vector`` :
{xrst_spell_off}
{xrst_code cpp} */
   typedef CppAD::sparse_rcv<s_vector, a1_vector> a1_sparse_rcv;
/* {xrst_code}
{xrst_spell_on}

End Namespace
*************
{xrst_spell_off}
{xrst_code cpp} */
} }
/* {xrst_code}
{xrst_spell_on}

{xrst_end typedef}
*/

# endif
