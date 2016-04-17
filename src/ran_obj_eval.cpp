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

/*
$begin ran_obj_eval$$
$spell
	chol
	hes
	CppAD
	ran_obj
	cppad
	obj
	eval
	vec
	const
	Cpp
	cholesky
$$

$section Evaluate Laplace Approximation and Random Objective$$

$head Syntax$$
$icode%h% = %mixed_object%.ran_obj_eval(%fixed_vec%, %random_vec%)%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

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

$head chol_ran_hes_$$
It is assumed that the member variable
$codei%
	CPPAD_MIXED_CHOLESKY chol_ran_hes_
%$$
was updated using $cref update_factor$$ for the specified values of the
fixed and random effects.

$head fixed_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fixed_vec%
%$$
It specifies the value of the
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$
vector $latex \theta$$ at which $latex h( \theta , u)$$ is evaluated.

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
It specifies the value of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$ at which $latex h( \theta , u)$$ is evaluated.
Note that the Laplace approximation is equal to the random objective when
$latex u$$ is the
$cref/optimal random effects
	/theory
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
	example/private/ran_obj_eval_xam.cpp
%$$
$head Example$$
The file $cref ran_obj_eval_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/

// ----------------------------------------------------------------------------
double cppad_mixed::ran_obj_eval(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec )
{	assert( init_ran_like_done_ );
	assert( init_ran_hes_done_ );

	assert( fixed_vec.size() == n_fixed_ );
	assert( random_vec.size() == n_random_ );

# ifndef NDEBUG
	// number of non-zeros in Hessian
	size_t K = ran_hes_.row.size();
	assert( K == ran_hes_.col.size() );
# endif

	// pack fixed and random effects into one vector
	d_vector both(n_fixed_ + n_random_);
	pack(fixed_vec, random_vec, both);

	// compute the logdet( f_{u,u}(theta, u )
	double logdet = chol_ran_hes_.logdet();

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


