// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_UNPACK_HPP
# define CPPAD_MIXED_UNPACK_HPP
# include <cppad/mixed/cppad_mixed.hpp>

/*
{xrst_begin unpack}

Pack Fixed Effect and Random Effects Into One Vector
####################################################

Syntax
******

| *mixed_object* . ``unpack`` ( *fixed_one* , *random_vec* , *both_vec* )
| *mixed_object* . ``unpack`` ( *fixed_one* , *fixed_two* , *random_vec* , *three_vec* )

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

Float_pack
**********
This can be any type.

Float_unpack
************
If *x* has type *Float_pack* ,
the syntax *Float_unpack* ( *x* ) must convert *x*
to the type *Float_unpack* .

fixed_one
*********
This argument has prototype

   ``CppAD::vector<`` *Float_unpack* >& *fixed_one*

The input value of its elements does not matter.
Upon return, it contains the value of the first
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>` .
The size of this vector must be equal to
:ref:`private_base_class@n_fixed_` .

fixed_two
*********
This argument has prototype

   ``CppAD::vector<`` *Float_unpack* >& *fixed_two*

The input value of its elements does not matter.
Upon return, it contains the value of the second
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>` .
The size of this vector must be equal to
:ref:`private_base_class@n_fixed_` .

random_vec
**********
This argument has prototype

   ``CppAD::vector<`` *Float_unpack* >& *random_vec*

The input value of its elements does not matter.
Upon return, it contains the value of the
:ref:`random effects<problem@Notation@Random Effects, u>` .
The size of this vector must be equal to
:ref:`private_base_class@n_random_` .

both_vec
********
This argument has prototype

   ``const CppAD::vector<`` *Float_pack* >& *both_vec*

If present, the size of this vector must be equal to
*n_fixed_* + *n_random_* .
It contains the fixed effect and random effects as one vector.

three_vec
*********
This argument has prototype

   ``const CppAD::vector<`` *Float_pack* >& *three_vec*

If present, the size of this vector must be equal to
2* *n_fixed_* + *n_random_* .
It contains the first fixed effects vector,
the second fixed effects vector,
and random effects as one vector.

{xrst_end unpack}
*/

template <class Float_unpack, class Float_pack>
void cppad_mixed::unpack(
   CppAD::vector<Float_unpack>&      fixed_one  ,
   CppAD::vector<Float_unpack>&      random_vec ,
   const CppAD::vector<Float_pack>&  both_vec   ) const
{
   assert( fixed_one.size() == n_fixed_ );
   assert( random_vec.size() == n_random_ );
   assert( both_vec.size() == n_fixed_ + n_random_ );
   for(size_t j = 0; j < n_fixed_; j++)
      fixed_one[j] = Float_unpack( both_vec[j] );
   for(size_t j = 0; j < n_random_; j++)
      random_vec[j] = Float_unpack( both_vec[n_fixed_ + j] );
}

template <class Float_unpack, class Float_pack>
void cppad_mixed::unpack(
   CppAD::vector<Float_unpack>&      fixed_one  ,
   CppAD::vector<Float_unpack>&      fixed_two  ,
   CppAD::vector<Float_unpack>&      random_vec ,
   const CppAD::vector<Float_pack>&  three_vec  ) const
{
   assert( fixed_one.size() == n_fixed_ );
   assert( fixed_two.size() == n_fixed_ );
   assert( random_vec.size() == n_random_ );
   assert( three_vec.size() == 2 * n_fixed_ + n_random_ );
   for(size_t j = 0; j < n_fixed_; j++)
   {  fixed_one[j] = Float_unpack( three_vec[j] );
      fixed_two[j] = Float_unpack( three_vec[n_fixed_ + j] );
   }
   for(size_t j = 0; j < n_random_; j++)
      random_vec[j] = Float_unpack( three_vec[2 * n_fixed_ + j] );
}



# endif
