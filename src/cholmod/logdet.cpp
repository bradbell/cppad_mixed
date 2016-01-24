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
$begin cholmod_logdet$$
$spell
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
The $cref cholmod$$ class is an
$cref/implementation detail/cholesky/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head cholmod_obj$$
This object has prototype
$codei%
	const CppAD::mixed::cholmod %cholmod_obj%
%$$
In addition, it must have a previous call to
$cref cholmod_update$$.

$head logdet$$
This return value has prototype
$codei%
	double %logdet%
%$$
Is the log of the determinant of the Hessian corresponding
to the previous call to $cref cholmod_update$$.

$head Example$$
The file $cref/cholmod_xam.cpp/cholmod_xam.cpp/logdet/$$ contains an
example and test that uses this function.

$end
*/

# include <cppad/mixed/cholmod.hpp>
# include <cassert>
# include <cmath>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

double cholmod::logdet(void) const
{	// factorization A = LDL'
	double logdet_A = 0.0;
	int*    L_p  = (int *) factor_->p;
	int*    L_i  = (int *) factor_->i;
	double* L_x  = (double *) factor_->x;
	for(size_t j = 0; j < nrow_; j++)
	{	// first element for each column is always the diagonal element
		assert( size_t( L_i [ L_p[j] ] ) == j );
		// j-th element on diagonal of D in factorization
		double dj = L_x[ L_p[j] ];
		assert( dj > 0.0 );
		logdet_A += std::log(dj);
	}
	return logdet_A;
}

} } // END_CPPAD_MIXED_NAMESPACE