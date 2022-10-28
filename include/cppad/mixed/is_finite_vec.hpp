// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_IS_FINITE_VEC_HPP
# define CPPAD_MIXED_IS_FINITE_VEC_HPP

/*
$begin is_finite_vec$$
$spell
	dev
	CppAD
	vec
$$

$section Are All Elements of a Vector Finite$$

$head Syntax$$
$icode%finite% = CppAD::mixed::is_finite_vec(%vec%)%$$

$head Prototype$$
$srcthisfile%
	0%// BEGIN_PROTOTYPE%// END_PROTOTYPE%1
%$$

$head Vector$$
Is a simple vector class.

$end
*/

// BEGIN_PROTOTYPE
namespace CppAD { namespace mixed {
	template <class Vector>
	bool is_finite_vec(const Vector& vec)
// END_PROTOTYPE
	{	typedef typename Vector::value_type scalar;
		scalar upper = scalar( std::numeric_limits<double>::infinity() );
		scalar lower = - upper;
		bool result  = true;
		size_t n     = vec.size();
		for(size_t i = 0; i < n; ++i)
		{	result &= lower < vec[i];
			result &= vec[i] < upper;
		}
		return result;
	}
} }

# endif
