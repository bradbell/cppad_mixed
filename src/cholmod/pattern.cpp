/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin ldlt_cholmod_pattern$$
$spell
	rcv
	rc
	ldlt_obj
	CppAD
	cholmod
	init
$$

$nospell
$bold This is cppad_mixed--20220519 documentation:$$ Here is a link to its
$href%https://cppad-mixed.readthedocs.io%current documentation%$$.
$$
$section Update Factorization Using new Matrix Values$$

$head Syntax$$
$icode%H_rc% = %ldlt_obj%.pattern()%$$

$head Prototype$$
$srcthisfile%0%// BEGIN_PROTOTYPE%// END_PROTOTYPE%1%$$

$head Private$$
The $cref ldlt_cholmod$$ class is an
$cref/implementation detail/ldlt_cholmod/Private/$$ and not part of the
CppAD Mixed user API.

$head ldlt_obj$$
This object has prototype
$codei%
	CppAD::mixed::ldlt_cholmod %ldlt_obj%
%$$
In addition, it must have a previous call to
$cref ldlt_cholmod_init$$.

$head H_rc$$
The return value is a copy of the sparsity pattern
$cref/H_rc/ldlt_cholmod_init/H_rc/$$ in the corresponding call to
$icode%ldlt_obj%.init(%H_rc%)%$$.

$head Order of Operations$$
This $icode ldlt_obj$$ function must be called,
after the constructor and $cref/init/ldlt_cholmod_init/$$.

$head Example$$
The file $cref/ldlt_cholmod.cpp/ldlt_cholmod.cpp/pattern/$$ contains an
example and test that uses this function.

$end
*/
// ----------------------------------------------------------------------------


# include <cppad/mixed/include_cholmod.hpp>
# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cassert>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN_PROTOTYPE
const sparse_rc& ldlt_cholmod::pattern(void) const
// END_PROTOTYPE
{	assert( init_done_ );
	return H_rc_;
}

} } // END_CPPAD_MIXED_NAMESPACE
