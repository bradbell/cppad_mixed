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
	cholmod_obj
	Cholesky
	logdet
	CppAD
$$

$section Compute Log Determinant for Current Cholesky Factor$$

$head Syntax$$
$icode%logdet% = %cholmod_obj%.logdet()
%$$

$head Private$$
The $cref ldlt_cholmod$$ class is an
$cref/implementation detail/ldlt_cholmod/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head cholmod_obj$$
This object has prototype
$codei%
	const CppAD::mixed::cholmod %cholmod_obj%
%$$
In addition, it must have a previous call to
$cref ldlt_cholmod_update$$.

$head logdet$$
This return value has prototype
$codei%
	double %logdet%
%$$
Is the log of the determinant of the Hessian corresponding
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

double cholmod::logdet(void) const
{	// factorization P H P' = L D L'
	double logdet_H = 0.0;
	int*    L_p  = (int *) factor_->p;
	double* L_x  = (double *) factor_->x;
# ifndef NDEBUG
	int*    L_i  = (int *) factor_->i;
# endif
	for(size_t j = 0; j < nrow_; j++)
	{	// first element for each column is always the diagonal element
		assert( size_t( L_i [ L_p[j] ] ) == j );
		// j-th element on diagonal of D in factorization
		double dj = L_x[ L_p[j] ];
		assert( dj > 0.0 );
		logdet_H += std::log(dj);
	}
	return logdet_H;
}

} } // END_CPPAD_MIXED_NAMESPACE
