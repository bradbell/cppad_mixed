// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/chol_hes_ran.hpp>

/*
$begin ranobj_eval$$
$spell
	CppAD
	ranobj
	cppad
	obj
	eval
	vec
	const
	Cpp
$$

$section Evaluate Laplace Approximation and Random Objective$$

$head Syntax$$
$icode%h% = %mixed_object%.ranobj_eval(%fixed_vec%, %random_vec%)%$$

$head Purpose$$
This routine evaluates the Laplace approximation
$cref/h(theta, u)
	/theory
	/Objective
	/Laplace Approximation, h(theta, u)
/$$.
Note that if the random effects are optimal,
then the Laplace approximation is equal to the random objective.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head fixed_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fixed_vec%
%$$
It specifies the value of the
$cref/fixed effects/cppad_mixed/Fixed Effects, theta/$$
vector $latex \theta$$ at which $latex h( \theta , u)$$ is evaluated.

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
It specifies the value of the
$cref/random effects/cppad_mixed/Random Effects, u/$$
vector $latex u$$ at which $latex h( \theta , u)$$ is evaluated.
Note that the Laplace approximation is equal to the random objective when
$latex u$$ is the
$cref/optimal random effects
	/theory
	/Objective
	/Optimal Random Effects, u^(theta)
/$$
$latex \hat{u} ( \theta )$$.

$head h$$
The return value has prototype
$codei%
	double %h%
%$$
and is the value of the Laplace approximation.

$children%
	example/private/ranobj_eval_xam.cpp
%$$
$head Example$$
The file $cref ranobj_eval_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/

// ----------------------------------------------------------------------------
double cppad_mixed::ranobj_eval(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec )
{	assert( init_ran_like_done_ );
	assert( init_hes_ran_done_ );

	assert( fixed_vec.size() == n_fixed_ );
	assert( random_vec.size() == n_random_ );

	// number of non-zeros in Hessian
	size_t K = hes_ran_.row.size();
	assert( K == hes_ran_.col.size() );

	// compute an LDL^T Cholesky factorization of f_{uu}^{(2)}(theta, u)
	d_vector both(n_fixed_ + n_random_);
	pack(fixed_vec, random_vec, both);
	CppAD::mixed::factorize_chol_hes_ran(
		n_fixed_, n_random_, hes_ran_.row, hes_ran_.col, both, hes_ran_fun_
	);

	// compute the logdet( f_{uu}^{(2)}(theta, u )
	double logdet = CppAD::mixed::logdet_chol_hes_ran(n_random_);

	// constant term
	double pi   = CppAD::atan(1.0) * 4.0;
	double constant_term = CppAD::log(2.0 * pi) * double(n_random_) / 2.0;

	// f(theta , u)
	d_vector vec = ran_like_fun_.Forward(0, both);
	assert( vec.size() == 1);

	// h(theta, u)
	double h = logdet / 2.0 + vec[0] - constant_term;

	return h;
}


