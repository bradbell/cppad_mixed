// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_UNPACK_HPP
# define CPPAD_MIXED_UNPACK_HPP
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/a2_double.hpp>

/*
$begin unpack$$
$spell
	CppAD
	cppad
	vec
	const
	Cpp
	dismod
	hpp
$$

$section Pack Fixed Effect and Random Effects Into One Vector$$

$head Syntax$$
$icode%mixed_object%.unpack(%fixed_one%, %random_vec%, %both_vec%)
%$$
$icode%mixed_object%.unpack(%fixed_one%, %fixed_two%, %random_vec%, %three_vec%)
%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head Float_pack$$
This can be any type.

$head Float_unpack$$
If $icode x$$ has type $icode Float_pack$$,
the syntax $icode%Float_unpack%(%x%)%$$ must convert $icode x$$
to the type $icode Float_unpack$$.

$head fixed_one$$
This argument has prototype
$codei%
	CppAD::vector<%Float_unpack%>& %fixed_one%
%$$
The input value of its elements does not matter.
Upon return, it contains the value of the first
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$.
The size of this vector must be equal to
$cref/n_fixed_/private/n_fixed_/$$.

$head fixed_two$$
This argument has prototype
$codei%
	CppAD::vector<%Float_unpack%>& %fixed_two%
%$$
The input value of its elements does not matter.
Upon return, it contains the value of the second
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$.
The size of this vector must be equal to
$cref/n_fixed_/private/n_fixed_/$$.

$head random_vec$$
This argument has prototype
$codei%
	CppAD::vector<%Float_unpack%>& %random_vec%
%$$
The input value of its elements does not matter.
Upon return, it contains the value of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$.
The size of this vector must be equal to
$cref/n_random_/private/n_random_/$$.

$head both_vec$$
This argument has prototype
$codei%
	const CppAD::vector<%Float_pack%>& %both_vec%
%$$
If present, the size of this vector must be equal to
$icode%n_fixed_% + %n_random_%$$.
It contains the fixed effect and random effects as one vector.

$head three_vec$$
This argument has prototype
$codei%
	const CppAD::vector<%Float_pack%>& %three_vec%
%$$
If present, the size of this vector must be equal to
$codei%2*%n_fixed_% + %n_random_%$$.
It contains the first fixed effects vector,
the second fixed effects vector,
and random effects as one vector.

$end
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
	{	fixed_one[j] = Float_unpack( three_vec[j] );
		fixed_two[j] = Float_unpack( three_vec[n_fixed_ + j] );
	}
	for(size_t j = 0; j < n_random_; j++)
		random_vec[j] = Float_unpack( three_vec[2 * n_fixed_ + j] );
}



# endif
