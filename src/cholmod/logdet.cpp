// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
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
$icode%logdet% = %ldlt_obj%.logdet(%sign%)
%$$

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

$head sign$$
This argument has prototype
$codei%
	int& %sign%
%$$
Its input value does no matter,
upon return it is $code +1$$ if the determinant is positive,
$code 0$$ if the determinant is zero,
and $code -1$$ if the determinant is negative.

$head logdet$$
This return value has prototype
$codei%
	double %logdet%
%$$
Is the log of the absolute value of the determinant corresponding
to the previous call to $cref ldlt_cholmod_update$$.

$head Example$$
The file $cref/ldlt_cholmod_xam.cpp/ldlt_cholmod_xam.cpp/logdet/$$ contains an
example and test that uses this function.

$end
*/

# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cassert>
# include <cmath>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

double ldlt_cholmod::logdet(int& sign) const
{	// factorization P H P' = L D L'
	double logdet_H = 0.0;
	int*    L_p  = (int *) factor_->p;
	double* L_x  = (double *) factor_->x;
# ifndef NDEBUG
	int*    L_i  = (int *) factor_->i;
# endif
	sign = 1;
	for(size_t j = 0; j < nrow_; j++)
	{	// first element for each column is always the diagonal element
		assert( size_t( L_i [ L_p[j] ] ) == j );
		// j-th element on diagonal of D in factorization
		double dj = L_x[ L_p[j] ];
		if( dj == 0.0 )
		{	sign = 0;
			return - std::numeric_limits<double>::infinity();
		}
		if( dj < 0.0 )
			sign = - sign;
		logdet_H += std::log( std::fabs(dj) );
	}
	return logdet_H;
}

} } // END_CPPAD_MIXED_NAMESPACE
