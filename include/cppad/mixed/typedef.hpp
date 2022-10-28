// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_TYPEDEF_HPP
# define CPPAD_MIXED_TYPEDEF_HPP

# include <cppad/cppad.hpp>

/*
$begin typedef$$
$spell
	CppAD
	Namespace
	typedef
	rc
	rcv
	cppad.hpp
$$

$section Types Defined in the CppAD Mixed Namespace$$

$head Syntax$$
$codei%# include <cppad/mixed/typedef.hpp>
%$$


$head Begin Namespace$$
All the definitions below are made inside the $code CppAD::mixed$$ namespace;
i.e.,
$srccode%cpp% */
namespace CppAD { namespace mixed {

/* %$$

$head Scalar Types$$

$subhead a1_double$$
Scalar with one level of AD:
$srccode%cpp% */
	typedef CppAD::AD<double> a1_double;
/* %$$

$head Vector Types$$

$subhead s_vector$$
Vectors with elements of type $code size_t$$:
$srccode%cpp% */
	typedef CppAD::vector<size_t> s_vector;
/* %$$

$subhead d_vector$$
Vectors with elements of type $code double$$:
$srccode%cpp% */
	typedef CppAD::vector<double> d_vector;
/* %$$

$subhead a1_vector$$
Vectors with elements of that have one level of AD:
$srccode%cpp% */
	typedef CppAD::vector<a1_double> a1_vector;
/* %$$

$head Sparse Types$$

$subhead sparse_rc$$
Sparsity patterns using index vector of type $code s_vector$$:
$srccode%cpp% */
	typedef CppAD::sparse_rc<s_vector> sparse_rc;
/* %$$

$subhead d_sparse_rcv$$
Sparse matrices using index vector of type $code s_vector$$
and value vectors of type $code d_vector$$:
$srccode%cpp% */
	typedef CppAD::sparse_rcv<s_vector, d_vector> d_sparse_rcv;

/* %$$
$subhead a1_sparse_rcv$$
Sparse matrices using index vector of type $code s_vector$$
and value vectors of type $code a1_vector$$:
$srccode%cpp% */
	typedef CppAD::sparse_rcv<s_vector, a1_vector> a1_sparse_rcv;
/* %$$

$head End Namespace$$
$srccode%cpp% */
} }
/* %$$
$end
*/

# endif
