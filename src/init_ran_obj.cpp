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
$begin init_ran_obj$$
$spell
	CppAD
	init
	ran_obj
	cppad
	obj
	vec
	const
	Cpp
$$

$section Second Order Representation of Random Objective$$

$head Syntax$$
$icode%mixed_object%.init_ran_obj(%fixed_vec%, %random_vec%)%$$

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

$head ran_obj_fun_$$
The input value of the member variable
$codei%
	CppAD::ADFun<double> ran_obj_fun_
%$$
does not matter.
Upon return it contains a second order accurate recording of the
approximate random objective; see
$cref/H(beta, theta, u)
	/theory
	/Approximate Random Objective, H(beta, theta, u)
/$$.

$end
*/
# include <Eigen/Sparse>
# include <cppad/mixed/cppad_mixed.hpp>



// ----------------------------------------------------------------------------
void cppad_mixed::init_ran_obj(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec )
{	assert( ! init_ran_obj_done_ );
	assert( record_newton_atom_done_ );

	//	create an a1d_vector containing (beta, theta, u)
	a1d_vector beta_theta_u( 2 * n_fixed_ + n_random_ );
	pack(fixed_vec, fixed_vec, random_vec, beta_theta_u);

	// start recording a1_double operations
	CppAD::Independent( beta_theta_u );

	// split back out to beta, theta, u
	a1d_vector beta(n_fixed_), theta(n_fixed_), u(n_random_);
	unpack(beta, theta, u, beta_theta_u);

	// evaluate gradient f_u (beta , u )
	a1d_vector grad(n_random_);
	grad = ran_like_grad(beta, u);

	// Evaluate the log determinant of f_{u,u} ( theta , u)
	// and Newton step s = f_{u,u} ( theta , u) f_u (beta, u)
	a1d_vector theta_u_v(n_fixed_ + 2 * n_random_ );
	for(size_t j = 0; j < n_fixed_; j++)
		theta_u_v[j] = theta[j];
	for(size_t j = 0; j < n_random_; j++)
	{	theta_u_v[n_fixed_ + j]             = u[j];
		theta_u_v[n_fixed_ + n_random_ + j] = grad[j];
	}
	a1d_vector logdet_step(1 + n_random_);
	newton_atom_.eval(theta_u_v, logdet_step);
	//
	// constant term
	double pi   = CppAD::atan(1.0) * 4.0;
	double constant_term = CppAD::log(2.0 * pi) * double(n_random_) / 2.0;
	//
	a1d_vector both(n_fixed_ + n_random_), f(1), H(1);
	// -----------------------------------------------------------------------
	// U(beta, theta, u)
	a1d_vector U(n_random_);
	for(size_t j = 0; j < n_random_; j++)
		U[j] = u[j] - logdet_step[1 + j];

	// evaluate gradient f_u (beta , U )
	grad = ran_like_grad(beta, U);

	// Evaluate the log determinant and newton step
	a1d_vector beta_U_v(n_fixed_ + 2 * n_random_ );
	for(size_t j = 0; j < n_fixed_; j++)
		beta_U_v[j] = beta[j];
	for(size_t j = 0; j < n_random_; j++)
	{	beta_U_v[n_fixed_ + j]             = U[j];
		beta_U_v[n_fixed_ + n_random_ + j] = grad[j];
	}
	newton_atom_.eval(beta_U_v, logdet_step);
	// -----------------------------------------------------------------------
	// W(beta, theta, u)
	a1d_vector W(n_random_);
	for(size_t j = 0; j < n_random_; j++)
		W[j] = U[j] - logdet_step[1 + j];

	// Evaluate the log determinant
	a1d_vector beta_W_v(n_fixed_ + 2 * n_random_ );
	for(size_t j = 0; j < n_fixed_; j++)
		beta_W_v[j] = beta[j];
	for(size_t j = 0; j < n_random_; j++)
	{	beta_W_v[n_fixed_ + j]             = W[j];
		beta_W_v[n_fixed_ + n_random_ + j] = 0.0;
	}
	newton_atom_.eval(beta_W_v, logdet_step);
	pack(beta, W, both);
	f    = ran_likelihood_a1fun_.Forward(0, both);
	H[0] = logdet_step[0] / 2.0 + f[0] - constant_term;
	//
	ran_obj_fun_.Dependent(beta_theta_u, H);
# ifdef NDEBUG
	ran_obj_fun_.optimize();
# endif
	//
	init_ran_obj_done_ = true;
	return;
}

