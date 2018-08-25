// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin ldlt_cholmod_logdet$$
$spell
	ldlt
	xam
	const
	ldlt_obj
	cholmod
	logdet
	CppAD
$$

$section Compute Log Determinant for Current Factor$$

$head Syntax$$
$icode%logdet% = %ldlt_obj%.logdet(%negative%)
%$$

$head Prototype$$
$srcfile%src/cholmod/logdet.cpp
	%0%// BEGIN_PROTOTYPE%// END_PROTOTYPE%1%$$

$head Private$$
The $cref ldlt_cholmod$$ class is an
$cref/implementation detail/ldlt_cholmod/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head ldlt_obj$$
This object has prototype
$codei%
	const CppAD::mixed::ldlt_cholmod %ldlt_obj%
%$$
In addition, it must have a previous call to
$cref ldlt_cholmod_update$$.

$head negative$$
The input value of this argument does no matter.
Upon return, it is the number of elements of
$cref/D/ldlt_cholmod/Factorization/D/$$
that are less than zero.

$head logdet$$
This return value
is the log of the absolute value of the determinant corresponding
to the previous call to $cref ldlt_cholmod_update$$.
If the matrix is singular, $icode logdet$$ is
minus infinity.

$head Example$$
The file $cref/ldlt_cholmod.cpp/ldlt_cholmod.cpp/logdet/$$ contains an
example and test that uses this function.

$end
*/

# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cassert>
# include <cmath>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN_PROTOTYPE
double ldlt_cholmod::logdet(size_t& negative) const
// END_PROTOTYPE
{	assert( update_called_ );
	//
	// factorization P H P' = L D L'
	int*    L_p  = (int *) factor_->p;
	double* L_x  = (double *) factor_->x;
# ifndef NDEBUG
	int*    L_i  = (int *) factor_->i;
# endif
	negative        = 0;
	bool has_zero   = false;
	double logdet_H = 0.0;
	for(size_t j = 0; j < nrow_; j++)
	{	// first element for each column is always the diagonal element
		assert( size_t( L_i [ L_p[j] ] ) == j );
		// j-th element on diagonal of D in factorization
		double dj = L_x[ L_p[j] ];
		has_zero |= dj == 0.0;
		if( dj < 0.0 )
			negative++;
		if( ! has_zero )
			logdet_H += std::log( std::fabs(dj) );
	}
	if( has_zero )
		return - std::numeric_limits<double>::infinity();
	//
	return logdet_H;
}

} } // END_CPPAD_MIXED_NAMESPACE
