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
	//
	assert( A_rcv_.nnz() == A_rcv_.row().size() );
	assert( A_rcv_.nnz() == A_rcv_.col().size() );
	assert( A_rcv_.nnz() == A_rcv_.val().size() );
	//
	//	beta
	a1_vector beta(n_fixed_);
	for(size_t j = 0; j < n_fixed_; ++j)
		beta[j] = fixed_vec[j];
	//
	// theta_u
	a1_vector theta_u( n_fixed_ + n_random_ );
	pack(fixed_vec, random_vec, theta_u);
	//
	// theta, u
	a1_vector theta(n_fixed_), u(n_random_);
	unpack(theta, u, theta_u);
	//
	// sparsity pattern with respect to (theta, u)
	// for the lower triangle of the random Hessian
	const sparse_rc& ran_hes_both_rc( ran_hes_rcv_.pat() );
	//
	// nnz
	size_t nnz = ran_hes_both_rc.nnz();
	//
	// ran_hes_uu_rc
	// ran_hes_both_rc
	sparse_rc ran_hes_uu_rc( n_random_, n_random_,            nnz);
	sparse_rc ran_hes_mix_rc(n_random_, n_fixed_ + n_random_, nnz);
	for(size_t k = 0; k < nnz; k++)
	{	size_t r = ran_hes_both_rc.row()[k];
		size_t c = ran_hes_both_rc.col()[k];
		//
		assert( r >= n_fixed_ );
		assert( c >= n_fixed_ );
		assert( c <= r );
		//
		// row and column relative to just random effect
		ran_hes_uu_rc.set(k,  r - n_fixed_, c - n_fixed_ );
		//
		// row relative to random effects, column relative to both
		ran_hes_mix_rc.set(k, r - n_fixed_, c);
	}
	// row, col
	const s_vector& row( ran_hes_uu_rc.row() );
	const s_vector& col( ran_hes_uu_rc.col() );
	//
	// hes_val
	a1_vector hes_val = ran_likelihood_hes(theta, u, row, col);
	if( hes_val.size() == 0 )
	{	// The user has not defined ran_likelihood_hes, so use AD to calcuate
		// the Hessian of the random likelihood w.r.t the random effects.
		hes_val.resize(nnz);
		a1_sparse_rcv subset( ran_hes_mix_rc );
		ran_jac_a1fun_.subgraph_jac_rev(theta_u, subset);
		ran_jac_a1fun_.clear_subgraph();
		hes_val = subset.val();
	}
	assert( hes_val.size() == nnz );
	//
	// ran_hes_uu_rcv
	a1_sparse_rcv ran_hes_uu_rcv( ran_hes_uu_rc );
	for(size_t k = 0; k < nnz; ++k)
		ran_hes_uu_rcv.set(k, hes_val[k] );
	//
	// ldlt_obj
	CppAD::mixed::ldlt_eigen<a1_double> ldlt_obj(n_random_);
	ldlt_obj.init( ran_hes_uu_rcv.pat() );
	ldlt_obj.update( ran_hes_uu_rcv );
	//
	// start recording a1_double operations
	// beta:    independent variables
	// theta_u: dynamic parameters
	size_t abort_op_index = 0;
	bool record_compare   = false;
	CppAD::Independent(beta, abort_op_index, record_compare, theta_u);
	//
	// theta, u
	unpack(theta, u, theta_u);
	//
	// beta_theta_u
	a1_vector beta_theta_u(2 * n_fixed_ + n_random_);
	pack(beta, theta, u, beta_theta_u);
	//
	// W(beta, theta, u)
	a1_vector W = CppAD::mixed::order2random(
		*this,
		n_fixed_,
		n_random_,
		ran_jac_a1fun_,
		ran_hes_both_rc,
		beta_theta_u
	);
	// -----------------------------------------------------------------------
	// Evaluate f_{uu} (beta , W).
	//
	// beta_W
	a1_vector beta_W(n_fixed_ + n_random_);
	pack(beta, W, beta_W);
	//
	// hes_val
	hes_val.resize(0); // so assignment works for both cases
	hes_val = ran_likelihood_hes(beta, W, row, col);
	if( hes_val.size() == 0 )
	{	// The user has not defined ran_likelihood_hes, so use AD to calcuate
		// the Hessian of the random likelihood w.r.t the random effects.
		hes_val.resize(nnz);
		a1_sparse_rcv subset( ran_hes_mix_rc );
		ran_jac_a1fun_.subgraph_jac_rev(beta_W, subset);
		ran_jac_a1fun_.clear_subgraph();
		hes_val = subset.val();
	}
	assert( hes_val.size() == nnz );
	//
	// ran_hes_uu_rcv.val() = hes_val
	for(size_t k = 0; k < nnz; ++k)
		ran_hes_uu_rcv.set(k, hes_val[k]);
	//
	// ldlt_obj.update
	ldlt_obj.update( ran_hes_uu_rcv );
	//
	// logdet [ f_{uu} ( beta , W ) ]
	size_t negative;
	a1_double logdet = ldlt_obj.logdet(negative);
	assert( negative == 0 );
	//
	// Evaluate the random likelihood using (beta, W)
	a1_vector f(1);
	f  = ran_like_a1fun_.Forward(0, beta_W);
	if( CppAD::hasnan(f) ) throw CppAD::mixed::exception(
		"init_laplace_obj", "result has a nan"
	);
	//
	// constant term
	double pi   = CppAD::atan(1.0) * 4.0;
	double constant_term = CppAD::log(2.0 * pi) * double(n_random_) / 2.0;
	//
	// now the random part of the Laplace objective
	a1_vector HB(1 + A_rcv_.nr());
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
	laplace_obj_fun_.Dependent(beta, HB);
	laplace_obj_fun_.check_for_nan(false);
# if CPPAD_MIXED_OPTIMIZE_CPPAD_FUNCTION
	laplace_obj_fun_.optimize("no_conditional_skip");
# endif
	//
	init_laplace_obj_done_ = true;
	return;
}
