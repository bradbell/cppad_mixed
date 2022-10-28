// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
$begin order2random$$
$spell
	CppAD
	jac
	cppad
	hes
	Jacobian
	ldlt
	uu
	const
	rc
$$

$section Second Order Representation of Random Effects, W(beta, theta, u)$$

$head Syntax$$
$icode%W% = CppAD::mixed::order2random(
	%n_fixed%,
	%n_random%,
	%jac_a1fun%,
	%a1_ldlt_ran_hes%,
	%beta%,
	%theta_u%
)%$$

$head Prototype$$
$srcthisfile%0%// BEGIN PROTOTYPE%// END PROTOTYPE%1%$$

$head Private$$
This $code cppad_mixed$$ is a $cref private_base_class$$ member function.

$head n_fixed$$
number of fixed effects.

$head n_random$$
number of random effects

$head jac_a1fun$$
This is the Jacobian of the random likelihood
with respect to the random effects
$latex f_u ( \theta, u )$$.
The domain indices (range indices) are with respect to the random effects
(fixed and random effects); i.e.,
$codei%
	%jac_a1fun%.Domain() == %n_fixed% + %n_random%
	%jac_a1fun%.Range() == %n_random%
%$$

$head a1_ldlt_ran_hes$$
This is the LDLT factorization for the
for the Hessian w.r.t the random effect; i.e., f_uu ( theta , u ).
Note that this object is $code const$$ and the previous call to
$cref/update/ldlt_eigen_update/$$ corresponded to (theta, u).

$head beta$$
This vector has size $icode n_fixed$$ and is the CppAD independent
variable vector at which we are evaluating the second order approximation
for the optimal random effects
$cref/W(beta, theta, u)
	/theory
	/Approximate Optimal Random Effects
	/Second Order, W(beta, theta, u)
/$$.

$head theta_u$$
This vector has size $icode%n_fixed% + %n_random%$$
and is the CppAD independent dynamic parameter vector
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

namespace { // BEGIN_EMPTY_NAMESPACE

typedef CppAD::mixed::a1_double                     a1_double;
typedef Eigen::Matrix<a1_double, Eigen::Dynamic, 1> a1_eigen_vector;

} // END_EMPTY_NAMESPACE


namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN PROTOTYPE
a1_vector order2random(
	size_t                              n_fixed         ,
	size_t                              n_random        ,
	CppAD::ADFun<a1_double, double>&    jac_a1fun       ,
	const ldlt_eigen<a1_double>&        a1_ldlt_ran_hes ,
	const a1_vector&                    beta            ,
	const a1_vector&                    theta_u         )
// END PROTOTYPE
{	assert( beta.size() == n_fixed );
	assert( theta_u.size() == n_fixed + n_random );
	assert( jac_a1fun.Domain() == n_fixed + n_random );
	assert( jac_a1fun.Range()  == n_random);
	//
	// L * D * L^T = P * H * P^T
	Eigen::SparseMatrix<a1_double, Eigen::ColMajor> L;
	a1_eigen_vector                                 D;
	Eigen::PermutationMatrix<Eigen::Dynamic>        P;
	a1_ldlt_ran_hes.split(L, D, P);
	// b  = H * x
	// b  = P^T * L * D * L^T * P * x
	// -----------------------------------------------------------------------
	// first partial Newton step
	//------------------------------------------------------------------------
	//
	// beta_u
	a1_vector beta_u(n_fixed + n_random);
	for(size_t j = 0; j < n_fixed; ++j)
		beta_u[j] = beta[j];
	for(size_t j = 0; j < n_random; ++j)
		beta_u[j + n_fixed] = theta_u[j + n_fixed];
	//
	// grad = f_u (beta , u )
	a1_vector grad(n_random);
	grad = jac_a1fun.Forward(0, beta_u);
	//
	// grad = f_{u,u} (theta, u) * step
	a1_eigen_vector step(n_random);
	for(size_t j = 0; j < n_random; ++j)
		step[j] = grad[j];
	step = ldlt_eigen<a1_double>::solve_LDLT(L, D, P, step);
	//
	// U(beta, theta, u) = u - step
	a1_vector U(n_random);
	for(size_t j = 0; j < n_random; j++)
		U[j] = theta_u[j + n_fixed] - step[j];
	// -----------------------------------------------------------------------
	// second partial Newton step
	//------------------------------------------------------------------------
	// grad = f_u (beta , U )
	for(size_t j = 0; j < n_random; ++j)
		beta_u[n_fixed + j] = U[j];
	grad = jac_a1fun.Forward(0, beta_u);
	//
	// grad = f_{u,u} (theta, u) * step
	for(size_t j = 0; j < n_random; ++j)
		step[j] = grad[j];
	step = ldlt_eigen<a1_double>::solve_LDLT(L, D, P, step);
	//
	// W(beta, theta, u) = U - step
	a1_vector W(n_random);
	for(size_t j = 0; j < n_random; j++)
		W[j] = U[j] - step[j];
	//------------------------------------------------------------------------
	return W;
}

} } // END_CPPAD_MIXED_NAMESPACE
