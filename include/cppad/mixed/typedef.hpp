// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_TYPEDEF_HPP
# define CPPAD_MIXED_TYPEDEF_HPP

# include <cppad/cppad.hpp>

/*
$begin cppad_mixed_typedef$$
$spell
	CppAD
	Namespace
	typedef
	rcv
	cppad.hpp
$$

$section Types Defined in the CppAD Mixed Namespace$$

$head Syntax$$
$codei%# include <cppad/mixed/typedef.hpp>
%$$
Note that this is not necessary if you include
$code cppad/mixed/cppad_mixed.hpp$$.


$head Scalar Types$$

$subhead a1_double$$
Scalar with one level of AD:
$srccode%cpp% */
namespace CppAD { namespace mixed {
	typedef CppAD::AD<double> a1_double;
} }
/* %$$

$subhead a2_double$$
Scalar with two levels of AD:
$srccode%cpp% */
namespace CppAD { namespace mixed {
	typedef CppAD::AD<a1_double> a2_double;
} }
/* %$$

$head Vector Types$$

$subhead s_vector$$
Vectors with elements of type $code size_t$$:
$srccode%cpp% */
namespace CppAD { namespace mixed {
	typedef CppAD::vector<size_t> s_vector;
} }
/* %$$

$subhead d_vector$$
Vectors with elements of type $code double$$:
$srccode%cpp% */
namespace CppAD { namespace mixed {
	typedef CppAD::vector<double> d_vector;
} }
/* %$$

$subhead a1_vector$$
Vectors with elements of that have one level of AD:
$srccode%cpp% */
namespace CppAD { namespace mixed {
	typedef CppAD::vector<a1_double> a1_vector;
} }
/* %$$

$subhead a2_vector$$
Vectors with elements of that have two levels of AD:
$srccode%cpp% */
namespace CppAD { namespace mixed {
	typedef CppAD::vector<a2_double> a2_vector;
} }
/* %$$

$head Sparse Matrices$$
Sparse matrices using index vector of type $code s_vector$$
and value vectors of type $code d_vector$$:
$srccode%cpp% */
namespace CppAD { namespace mixed {
	typedef CppAD::sparse_rcv<s_vector, d_vector> sparse_rcv;
} }
/* %$$
$end
*/

# endif
