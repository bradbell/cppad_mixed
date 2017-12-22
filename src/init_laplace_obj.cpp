// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-17 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin init_ran_objcon$$
$spell
	rcv
	nr
	objcon
	CppAD
	init
	ran_obj
	cppad
	obj
	vec
	const
	Cpp
$$

$section Second Order Representation of Laplace Objective and Constraints$$

$head Syntax$$
$icode%mixed_object%.init_ran_objcon(%fixed_vec%, %random_vec%)%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head Assumptions$$
The member variables
$cref/init_ran_like_done_/init_ran_like/init_ran_like_done_/$$ and
$cref/init_newton_checkpoint_done_/initialize/init_newton_checkpoint_done_/$$
are true.

$head init_laplace_obj_done_$$
The input value of this member variable must be false.
Upon return it is true.

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
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$
vector $latex \theta$$ at which the initialization is done.

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
It specifies the value of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$ at which the initialization is done.

$head laplace_obj_fun_$$
The input value of the member variable
$codei%
	CppAD::ADFun<double> laplace_obj_fun_
%$$
does not matter.
Upon return it contains a second order accurate recording of the
approximate Laplace objective; see
$cref/H(beta, theta, u)
	/theory
	/Approximate Laplace Objective, H(beta, theta, u)
/$$
$latex H( \beta , \theta , u )$$,
followed by the approximate random constraint function; see
$cref/B(beta, theta, u)
	/theory
	/Approximate Random Constraint Function, B(beta, theta, u)
/$$
$latex B( \beta , \theta , u )$$. Thus
$codei%
	laplace_obj_fun_.Domain() == n_fixed_ + n_fixed_ + n_random_
	laplace_obj_fun_.Range()  == 1 + A_rcv_.nr()
%$$


$end
*/
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/configure.hpp>

// ----------------------------------------------------------------------------
void cppad_mixed::init_ran_objcon(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec )
{	assert( ! init_laplace_obj_done_ );
	assert( init_ran_like_done_ );
	assert( init_newton_checkpoint_done_ );
	//
	assert( A_rcv_.nnz() == A_rcv_.row().size() );
	assert( A_rcv_.nnz() == A_rcv_.col().size() );
	assert( A_rcv_.nnz() == A_rcv_.val().size() );
	//
	// constant term
	double pi   = CppAD::atan(1.0) * 4.0;
	double constant_term = CppAD::log(2.0 * pi) * double(n_random_) / 2.0;
	//
	a1_vector f(1), HB(1 + A_rcv_.nr());
	//
	//	create an a1_vector containing (beta, theta, u)
	a1_vector beta_theta_u( 2 * n_fixed_ + n_random_ );
	pack(fixed_vec, fixed_vec, random_vec, beta_theta_u);
	//
	// start recording a1_double operations
	CppAD::Independent( beta_theta_u );
	//
	// split back out to beta, theta, u
	a1_vector beta(n_fixed_), theta(n_fixed_), u(n_random_);
	unpack(beta, theta, u, beta_theta_u);
	// -----------------------------------------------------------------------
	// First Newton Step
	//------------------------------------------------------------------------
	// evaluate gradient f_u (beta , u )
	a1_vector grad(n_random_);
	grad = ran_like_jac(beta, u);
	//
	// Evaluate the log determinant of f_{u,u} ( theta , u)
	// and Newton step s = f_{u,u} ( theta , u) f_u (beta, u)
	a1_vector theta_u_v(n_fixed_ + 2 * n_random_ );
	for(size_t j = 0; j < n_fixed_; j++)
		theta_u_v[j] = theta[j];
	for(size_t j = 0; j < n_random_; j++)
	{	theta_u_v[n_fixed_ + j]             = u[j];
		theta_u_v[n_fixed_ + n_random_ + j] = grad[j];
	}
	a1_vector logdet_step(1 + n_random_);
	newton_checkpoint_.eval(theta_u_v, logdet_step);
	//
	// U(beta, theta, u)
	a1_vector U(n_random_);
	for(size_t j = 0; j < n_random_; j++)
		U[j] = u[j] - logdet_step[1 + j];
	// -----------------------------------------------------------------------
	// Second Newton step
	//------------------------------------------------------------------------
	// evaluate gradient f_u (beta , U )
	grad = ran_like_jac(beta, U);
	//
	// Evaluate the log determinant and second newton step
	// s = f_{u,u} ( theta , u) f_u (beta, U)
	a1_vector theta_U_v(n_fixed_ + 2 * n_random_ );
	for(size_t j = 0; j < n_fixed_; j++)
		theta_U_v[j] = theta[j];
	for(size_t j = 0; j < n_random_; j++)
	{	theta_U_v[n_fixed_ + j]             = U[j];
		theta_U_v[n_fixed_ + n_random_ + j] = grad[j];
	}
	newton_checkpoint_.eval(theta_U_v, logdet_step);
	//
	// W(beta, theta, u)
	a1_vector W(n_random_);
	for(size_t j = 0; j < n_random_; j++)
		W[j] = U[j] - logdet_step[1 + j];
	// -----------------------------------------------------------------------
	//
	// Evaluate the log determinant using (beta, W)
	a1_vector beta_W_v(n_fixed_ + 2 * n_random_ );
	for(size_t j = 0; j < n_fixed_; j++)
		beta_W_v[j] = beta[j];
	for(size_t j = 0; j < n_random_; j++)
	{	beta_W_v[n_fixed_ + j]             = W[j];
		beta_W_v[n_fixed_ + n_random_ + j] = 0.0;
	}
	newton_checkpoint_.eval(beta_W_v, logdet_step);
	//
	// Evaluate the random likelihood using (beta, U)
	a1_vector beta_U(n_fixed_ + n_random_);
	pack(beta, U, beta_U);
	f     = ran_like_a1fun_.Forward(0, beta_U);
	if( CppAD::hasnan(f) ) throw CppAD::mixed::exception(
		"init_ran_objcon", "result has a nan"
	);
	//
	// now the random part of the Laplace objective
	HB[0] = logdet_step[0] / 2.0 + f[0] - constant_term;
	//
	if( A_rcv_.nr() > 0 )
	{
		// multiply the matrix A times the vector W and put result in
		// HB (with an index offset of 1)
		for(size_t i = 0; i < A_rcv_.nr(); i++)
			HB[1 + i] = a1_double(0.0);

		size_t K = A_rcv_.row().size();
		for(size_t k = 0; k < K; k++)
		{	size_t i = A_rcv_.row()[k];
			size_t j = A_rcv_.col()[k];
			double v = A_rcv_.val()[k];
			HB[1 + i]  += a1_double(v) * W[j];
		}
	}
	//
	laplace_obj_fun_.Dependent(beta_theta_u, HB);
	laplace_obj_fun_.check_for_nan(false);
# if CPPAD_MIXED_OPTIMIZE_CPPAD_FUNCTION
	laplace_obj_fun_.optimize("no_conditional_skip");
# endif
	//
	init_laplace_obj_done_ = true;
	return;
}
