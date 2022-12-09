// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_PACK_HPP
# define CPPAD_MIXED_PACK_HPP
# include <cppad/mixed/cppad_mixed.hpp>

/*
{xrst_begin pack dev}

Pack Fixed Effect and Random Effects Into One Vector
####################################################

Syntax
******

| *mixed_object* . ``pack`` ( *fixed_one* , *random_vec* , *both_vec* )
| *mixed_object* . ``pack`` ( *fixed_one* , *fixed_two* , *random_vec* , *three_vec* )

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

Float_unpack
************
This can be any type.

Float_pack
**********
If *x* has type *Float_unpack* ,
the syntax *Float_pack* ( *x* ) must convert *x*
to the type *Float_pack* .

fixed_one
*********
This argument has prototype

   ``const CppAD::vector<`` *Float_unpack* >& *fixed_one*

It specifies the a value for the
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>` .
The size of this vector must be equal to
:ref:`private_base_class@n_fixed_` .

fixed_two
*********
This argument has prototype

   ``const CppAD::vector<`` *Float_unpack* >& *fixed_two*

If present, it also specifies the a value for the
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>` .
The size of this vector must be equal to
:ref:`private_base_class@n_fixed_` .

random_vec
**********
This argument has prototype

   ``const CppAD::vector<`` *Float_unpack* >& *random_vec*

It specifies a value for the
:ref:`random effects<problem@Notation@Random Effects, u>` .
The size of this vector must be equal to
:ref:`n_fixed_<private_base_class@n_random_>` .

both_vec
********
This argument has prototype

   ``CppAD::vector<`` *Float_pack* >& *both_vec*

If present, the size of this vector must be equal to
*n_fixed_* + *n_random_* .
The input value of its elements does not matter.
Upon return, it contains the values in
*fixed_one* and *random_vec* as one vector in that order;
i.e., *fixed_one* comes first and then *random_vec* .

three_vec
*********
This argument has prototype

   ``CppAD::vector<`` *Float_pack* >& *three_vec*

If present, the size of this vector must be equal to
2* *n_fixed_* + *n_random_* .
The input value of its elements does not matter.
Upon return, it contains the values in
*fixed_one* , *fixed_two* , and *random_vec* as one vector.
The order of the result is unspecified.

{xrst_end pack}
*/

template <class Float_unpack, class Float_pack>
void cppad_mixed::pack(
   const CppAD::vector<Float_unpack>& fixed_one  ,
   const CppAD::vector<Float_unpack>& random_vec ,
   CppAD::vector<Float_pack>&         both_vec   ) const
{
   assert( fixed_one.size() == n_fixed_ );
   assert( random_vec.size() == n_random_ );
   assert( both_vec.size() == n_fixed_ + n_random_ );
   for(size_t j = 0; j < n_fixed_; j++)
      both_vec[j] = Float_pack( fixed_one[j] );
   for(size_t j = 0; j < n_random_; j++)
      both_vec[n_fixed_ + j] = Float_pack( random_vec[j] );
}
template <class Float_unpack, class Float_pack>
void cppad_mixed::pack(
   const CppAD::vector<Float_unpack>& fixed_one  ,
   const CppAD::vector<Float_unpack>& fixed_two  ,
   const CppAD::vector<Float_unpack>& random_vec ,
   CppAD::vector<Float_pack>&         three_vec  ) const
{
   assert( fixed_one.size() == n_fixed_ );
   assert( fixed_two.size() == n_fixed_ );
   assert( random_vec.size() == n_random_ );
   assert( three_vec.size() == 2 * n_fixed_ + n_random_ );
   for(size_t j = 0; j < n_fixed_; j++)
   {  three_vec[j] = Float_pack( fixed_one[j] );
      three_vec[n_fixed_ + j] = Float_pack( fixed_two[j] );
   }
   for(size_t j = 0; j < n_random_; j++)
      three_vec[2 * n_fixed_ + j] = Float_pack( random_vec[j] );
}



# endif
