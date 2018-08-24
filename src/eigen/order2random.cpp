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
$begin order2random$$
$spell
	CppAD
	jac
	cppad
	hes
	Jacobian
	rc
$$

$section Second Order Representation of Random Effects, W(beta, theta, u)$$

$head Syntax$$
$icode%W% = CppAD::mixed::order2random(
	%n_fixed%,
	%n_random%,
	%jac_a1fun%,
	%ran_hes_both_rc%,
	%beta_theta_u%
)%$$

$head Prototype$$
$srcfile%src/eigen/order2random.cpp
    %4%// BEGIN PROTOTYPE%// END PROTOTYPE%1%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head mixed_object$$
is the cppad_mixed object for this problem.
This object is only used to call
$codei%
	%mixed_object%.ran_likelihood_hes
%$$
and see if the user has defined the hessian of the random likelihood.

$head n_fixed$$
number of fixed effects.

$head n_random$$
number of random effects

$head jac_a1fun$$
This is the Jacobian of the random likelihood
with respect to the random effects
$latex f_u ( \theta, u )$$.

$head ran_hes_both_rc$$
This is the sparsity pattern for the lower triangle of the Hessian w.r.t the
random effect.
The indices in $icode ran_hes_both_rc$$ are with respect to both the fixed
and random effects and hence are greater than or equal $icode n_fixed$$.

$head beta_theta_u$$
This vector has size $codei%2*%n_fixed% + %n_random%$$
and specifies the value of $codei%(%beta%, %theta%, %u%)%$$
at which we are evaluating the second order approximation
for the optimal random effects
$cref/W(beta, theta, u)
	/theory
	/Approximate Optimal Random Effects
	/Second Order, W(beta, theta, u)
/$$.
$head W$$
The return value is $icode%W%(%beta%, %theta%, %u%)%$$.

$children%
	example/private/order2random.cpp
%$$
$head Example$$
The file $cref order2random.cpp$$ is an example and test of this
routine.

$end
*/

# include <Eigen/Sparse>
# include <cppad/mixed/order2random.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN PROTOTYPE
a1_vector order2random(
	cppad_mixed&             mixed_object    ,
	size_t                   n_fixed         ,
	size_t                   n_random        ,
	CppAD::ADFun<a1_double>& jac_a1fun       ,
	const sparse_rc&         ran_hes_both_rc ,
	const a1_vector&         beta_theta_u    )
// END PROTOTYPE
{	assert( beta_theta_u.size() == 2 * n_fixed + n_random );
	//
	// declare eigen matrix types
	typedef Eigen::Matrix<a1_double, Eigen::Dynamic, 1>     a1_eigen_vector;
	typedef Eigen::SparseMatrix<a1_double, Eigen::ColMajor> a1_eigen_sparse;
	//
	// beta, theta, u, theta_u, beta_u
	a1_vector beta(n_fixed), theta(n_fixed), u(n_random);
	a1_vector theta_u(n_fixed + n_random), beta_u(n_fixed + n_random);
	for(size_t j = 0; j < n_fixed; ++j)
	{	beta[j]    = beta_theta_u[j];
		beta_u[j]  = beta_theta_u[j];
		//
		theta[j]   = beta_theta_u[j + n_fixed];
		theta_u[j] = beta_theta_u[j + n_fixed];
	}
	for(size_t j = 0; j < n_random; ++j)
	{	u[j]                 = beta_theta_u[j + 2 * n_fixed ];
		theta_u[j + n_fixed] = beta_theta_u[j + 2 * n_fixed ];
		beta_u[j + n_fixed]  = beta_theta_u[j + 2 * n_fixed ];
	}
	// -----------------------------------------------------------------------
	// Evaluate f_{uu} (theta , u).
	//
	// n_low, row, col, val_out
	size_t n_low = ran_hes_both_rc.nnz();
	CppAD::vector<size_t> row(n_low), col(n_low);
	sparse_rc ran_hes_mix_rc(n_random, n_fixed + n_random, n_low);
	for(size_t k = 0; k < n_low; k++)
	{	assert( ran_hes_both_rc.row()[k] >= n_fixed );
		assert( ran_hes_both_rc.col()[k] >= n_fixed );
		//
		// row and column relative to just random effect
		row[k] = ran_hes_both_rc.row()[k] - n_fixed;
		col[k] = ran_hes_both_rc.col()[k] - n_fixed;
		//
		// row relative to random effects, coluimn relative to both
		ran_hes_mix_rc.set(k, row[k], col[k] + n_fixed);
	}
	a1_vector val_out = mixed_object.ran_likelihood_hes(theta, u, row, col);
	if( val_out.size() == 0 )
	{	// The user has not defined ran_likelihood_hes, so use AD to calcuate
		// the Hessian of the random likelihood w.r.t the random effects.
		val_out.resize(n_low);
		a1_sparse_rcv subset( ran_hes_mix_rc );
		jac_a1fun.subgraph_jac_rev(theta_u, subset);
		jac_a1fun.clear_subgraph();
		val_out = subset.val();
	}
	//
	// a1_hessian
	a1_eigen_sparse hessian;
	hessian.resize( int(n_random) , int(n_random) );
	for(size_t k = 0; k < n_low; k++)
		hessian.insert( int(row[k]), int(col[k]) ) = val_out[k];
	//
	// chol = L * D * L^T Cholesky factorization of f_{u,u} (theta, u)
	Eigen::SimplicialLDLT<a1_eigen_sparse, Eigen::Lower> chol;
	chol.analyzePattern(hessian);
	chol.factorize(hessian);
	// -----------------------------------------------------------------------
	// first partial Newton step
	//------------------------------------------------------------------------
	// grad = f_u (beta , u )
	a1_vector grad(n_random);
	grad = jac_a1fun.Forward(0, beta_u);
	//
	// step = f_{u,u} (theta, u)^{-1} * grad
	a1_eigen_vector step(n_random);
	for(size_t j = 0; j < n_random; ++j)
		step[j] = grad[j];
	step = chol.solve(step);
	//
	// U(beta, theta, u) = u - step
	a1_vector U(n_random);
	for(size_t j = 0; j < n_random; j++)
		U[j] = u[j] - step[j];
	// -----------------------------------------------------------------------
	// second partial Newton step
	//------------------------------------------------------------------------
	// grad = f_u (beta , U )
	for(size_t j = 0; j < n_random; ++j)
		beta_u[n_fixed + j] = U[j];
	grad = jac_a1fun.Forward(0, beta_u);
	//
	// step = f_{u,u} (theta, u)^{-1} * grad
	for(size_t j = 0; j < n_random; ++j)
		step[j] = grad[j];
	step = chol.solve(step);
	//
	// W(beta, theta, u) = U - step
	a1_vector W(n_random);
	for(size_t j = 0; j < n_random; j++)
		W[j] = U[j] - step[j];
	//------------------------------------------------------------------------
	return W;
}

} } // END_CPPAD_MIXED_NAMESPACE
