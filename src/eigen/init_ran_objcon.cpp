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
$begin init_ran_objcon$$
$spell
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

$section Second Order Representation of Random Objective and Constraints$$

$head Syntax$$
$icode%mixed_object%.init_ran_objcon(%fixed_vec%, %random_vec%)%$$

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

$head ran_objcon_fun_$$
The input value of the member variable
$codei%
	CppAD::ADFun<double> ran_objcon_fun_
%$$
does not matter.
Upon return it contains a second order accurate recording of the
approximate random objective; see
$cref/H(beta, theta, u)
	/theory
	/Approximate Random Objective, H(beta, theta, u)
/$$
$latex H( \beta , \theta , u )$$,
followed by the approximate random constraint function; see
$cref/B(beta, theta, u)
	/theory
	/Approximate Random Constraint Function, B(beta, theta, u)
/$$
$latex B( \beta , \theta , u )$$. Thus
$codei%
	ran_objcon_fun_.Domain() == n_fixed_ + n_fixed_ + n_random_
	ran_objcon_fun_.Range()  == 1 + n_ran_con_
%$$


$end
*/
# include <Eigen/Sparse>
# include <cppad/mixed/cppad_mixed.hpp>

// ----------------------------------------------------------------------------
void cppad_mixed::init_ran_objcon(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec )
{	assert( init_ran_con_done_ );
	assert( init_newton_atom_done_ );
	assert( ! init_ran_objcon_done_ );

	// number of constraints
	assert( n_ran_con_ == size_t ( ran_con_mat_.rows() ) );

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
	grad = ran_like_jac(beta, u);

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
	a1d_vector both(n_fixed_ + n_random_), f(1), HB(1 + n_ran_con_);
	// -----------------------------------------------------------------------
	// U(beta, theta, u)
	a1d_vector U(n_random_);
	for(size_t j = 0; j < n_random_; j++)
		U[j] = u[j] - logdet_step[1 + j];

	// evaluate gradient f_u (beta , U )
	grad = ran_like_jac(beta, U);

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

	// Evaluate the log determinant and random part of objective
	a1d_vector beta_W_v(n_fixed_ + 2 * n_random_ );
	for(size_t j = 0; j < n_fixed_; j++)
		beta_W_v[j] = beta[j];
	for(size_t j = 0; j < n_random_; j++)
	{	beta_W_v[n_fixed_ + j]             = W[j];
		beta_W_v[n_fixed_ + n_random_ + j] = 0.0;
	}
	newton_atom_.eval(beta_W_v, logdet_step);
	pack(beta, W, both);
	f     = ran_like_a1fun_.Forward(0, both);
	HB[0] = logdet_step[0] / 2.0 + f[0] - constant_term;

	if( n_ran_con_ > 0 )
	{
		// Copy W to an a1 eigen matrix
		using Eigen::Dynamic;
		typedef Eigen::Matrix<a1_double, Dynamic, Dynamic>  a1_eigen_dense;
		a1_eigen_dense eigen_W(n_random_, 1);
		for(size_t j = 0; j < n_random_; j++)
			eigen_W(j, 0) = W[j];

		// Copy ran_con_mat_ to an a1 eigen matrix
		using Eigen::ColMajor;
		typedef Eigen::SparseMatrix<double,    ColMajor>   eigen_sparse;
		typedef Eigen::SparseMatrix<a1_double, ColMajor>   a1_eigen_sparse;
		typedef typename eigen_sparse::InnerIterator       column_itr;
		a1_eigen_sparse eigen_A(n_ran_con_, n_random_);
		for(size_t j = 0; j < n_fixed_; j++)
		{	for(column_itr itr(ran_con_mat_, j); itr; ++itr)
			{	size_t i = itr.row();
				eigen_A.insert(i, j) = a1_double( itr.value() );
			}
		}

		// do the multiplication
		a1_eigen_dense eigen_B = eigen_A * eigen_W;

		// copy results to range vector
		for(size_t i = 0; i < n_ran_con_; i++)
			HB[1 + i] = eigen_B(i, 0);
	}
	//
	ran_objcon_fun_.Dependent(beta_theta_u, HB);
# ifdef NDEBUG
	ran_objcon_fun_.optimize();
# endif
	//
	init_ran_objcon_done_ = true;
	return;
}