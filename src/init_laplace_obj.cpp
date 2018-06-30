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
$begin init_laplace_obj$$
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
$icode%mixed_object%.init_laplace_obj(%fixed_vec%, %random_vec%)%$$

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
# include <Eigen/Sparse>
# include <cppad/mixed/order2random.hpp>
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/configure.hpp>

// ----------------------------------------------------------------------------
void cppad_mixed::init_laplace_obj(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec )
{	assert( ! init_laplace_obj_done_ );
	assert( init_ran_like_done_ );
	assert( init_newton_checkpoint_done_ );
	//
	// declare eigen matrix types
	typedef Eigen::Matrix<a1_double, Eigen::Dynamic, 1>     a1_eigen_vector;
	typedef Eigen::SparseMatrix<a1_double, Eigen::ColMajor> a1_eigen_sparse;
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
	size_t abort_op_index = 0;
	bool record_compare   = false;
	CppAD::Independent(beta_theta_u, abort_op_index, record_compare);
	//
	// split out beta, theta, u
	a1_vector beta(n_fixed_), theta(n_fixed_), u(n_random_);
	for(size_t j = 0; j < n_fixed_; ++j)
	{	beta[j]    = beta_theta_u[j];
		theta[j]   = beta_theta_u[j + n_fixed_];
	}
	for(size_t j = 0; j < n_random_; ++j)
		u[j]       = beta_theta_u[j + 2 * n_fixed_ ];
	//
	// W(beta, theta, u)
	const sparse_rc& ran_hes_rc( ran_hes_rcv_.pat() );
	a1_vector W = CppAD::mixed::order2random(
		*this,
		n_fixed_,
		n_random_,
		ran_like_a1fun_,
		ran_jac_a1fun_,
		ran_hes_rc,
		beta_theta_u
	);
	// -----------------------------------------------------------------------
	// beta_W
	a1_vector beta_W(n_fixed_ + n_random_);
	for(size_t j = 0; j < n_fixed_; j++)
		beta_W[j] = beta[j];
	for(size_t j = 0; j < n_random_; j++)
		beta_W[j + n_fixed_] = W[j];
	// -----------------------------------------------------------------------
	// Evaluate f_{uu} (beta , W).
	//
	// n_low, row, col, val_out
	size_t n_low = ran_hes_rc.nnz();
	CppAD::vector<size_t> row(n_low), col(n_low);
	sparse_rc ran_hes_mix_rc(n_random_, n_fixed_ + n_random_, n_low);
	for(size_t k = 0; k < n_low; k++)
	{	assert( ran_hes_rc.row()[k] >= n_fixed_ );
		assert( ran_hes_rc.col()[k] >= n_fixed_ );
		//
		// row and column relative to just random effect
		row[k] = ran_hes_rc.row()[k] - n_fixed_;
		col[k] = ran_hes_rc.col()[k] - n_fixed_;
		//
		// row relative to random effects, coluimn relative to both
		ran_hes_mix_rc.set(k, row[k], col[k] + n_fixed_);
	}
	a1_vector val_out = ran_likelihood_hes(beta, W, row, col);
	if( val_out.size() == 0 )
	{	// The user has not defined ran_likelihood_hes, so use AD to calcuate
		// the Hessian of the random likelihood w.r.t the random effects.
		val_out.resize(n_low);
		a1_sparse_rcv subset( ran_hes_mix_rc );
		ran_jac_a1fun_.subgraph_jac_rev(beta_W, subset);
		ran_jac_a1fun_.clear_subgraph();
		val_out = subset.val();
	}
	//
	// a1_hessian
	a1_eigen_sparse hessian;
	hessian.resize( int(n_random_) , int(n_random_) );
	for(size_t k = 0; k < n_low; k++)
		hessian.insert( int(row[k]), int(col[k]) ) = val_out[k];
	//
	// chol = L * D * L^T Cholesky factorization of f_{u,u} (beta, W)
	Eigen::SimplicialLDLT<a1_eigen_sparse, Eigen::Lower> chol;
	chol.analyzePattern(hessian);
	chol.factorize(hessian);
	// -----------------------------------------------------------------------
	//
	// logdet [ f_{uu} ( beta , W ) ]
	a1_double logdet     = 0.0;
	a1_eigen_vector diag = chol.vectorD();
	for(size_t j = 0; j < n_random_; ++j)
		logdet += log( diag[j] );
	//
	// Evaluate the random likelihood using (beta, W)
	f  = ran_like_a1fun_.Forward(0, beta_W);
	if( CppAD::hasnan(f) ) throw CppAD::mixed::exception(
		"init_laplace_obj", "result has a nan"
	);
	//
	// now the random part of the Laplace objective
	HB[0] = logdet / 2.0 + f[0] - constant_term;
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
