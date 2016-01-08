// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/chol_hes_ran.hpp>
/*
$begin update_factor$$
$spell
	Taylor
	vec
	Cholesky
	cppad
	const
	CppAD
	chol
	hes
$$

$section Update the Factorization of Hessian w.r.t. Random Effects$$

$head Syntax$$
$icode%mixed_object%.update_factor(%fixed_vec%, %random_vec%)%$$

$head Purpose$$
This routine updates the Cholesky factorization of the Hessian
$latex f_{u,u} ( \theta , u )^{-1} $$
so that it corresponds to the current value of
$latex \theta$$ and $latex u$$.

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head fixed_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fixed_vec%
%$$
and is the value of fixed effects $latex \theta$$.

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
and is the value of fixed effects $latex u$$.
If this factorization is used for computing derivatives of the
$cref/random objective/theory/Objective/Random Objective, r(theta)/$$,
this should be the
$cref/optimal random effects/theory/Optimal Random Effects, u^(theta)/$$.

$head hes_ran_fun_$$
The $cref/hes_ran_fun_/private/hes_ran_fun_/$$ member variable
will hold the first order Taylor coefficient corresponding
to the specified fixed and random effects; i.e.,
$codei%
	hes_ran_fun_.Forward(0, %both%)
%$$
has been called where $icode both$$ is a packed version
of the fixed and random effects.


$head factorize_chol_hes_ran_$$
On input, the static variable
$codei%
	CppAD::mixed::chol_hes_ran_
%$$
has been
$cref/analyzed/chol_hes_ran/analyze_chol_hes_ran/$$
using the sparsity pattern for the Hessian.
Upon return, $code chol_hes_ran_$$ contains the
$cref/factorization/chol_hes_ran/factorize_chol_hes_ran/$$
corresponding to the specified values for the fixed
and random effects.

$children%
	example/private/update_factor_xam.cpp
%$$
$head Example$$
The file $cref update_factor_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/
// ----------------------------------------------------------------------------
void cppad_mixed::update_factor(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec )
{	assert( init_hes_ran_done_ );
	assert( fixed_vec.size() == n_fixed_ );
	assert( random_vec.size() == n_random_ );

	// compute an LDL^T Cholesky factorization of f_{u,u} (theta, u)
	d_vector both(n_fixed_ + n_random_);
	pack(fixed_vec, random_vec, both);
	CppAD::mixed::factorize_chol_hes_ran(
		n_fixed_, n_random_, hes_ran_.row, hes_ran_.col, both, hes_ran_fun_
	);
}

