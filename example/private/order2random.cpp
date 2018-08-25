// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin order2random.cpp$$
$spell
	CppAD
	cppad
	hes
	interp
	xam
$$

$section order2random: Example and Test$$


$head Private$$
This example is not part of the
$cref/cppad_mixed public API/public/$$.

$code
$srcfile%example/private/order2random.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/mixed/order2random.hpp>

bool order2random_xam(void)
{
	bool   ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();
	//
	using CppAD::vector;
	typedef CppAD::AD<double>    a1_double;
	typedef CppAD::AD<a1_double> a2_double;
	//
	// n_fixed
	size_t n_fixed     = 3;
	//
	// n_random
	size_t n_random    = 3;
	//
	// n_both
	size_t n_both      = n_fixed + n_random;
	//
	// n_all
	size_t n_all       = 2 * n_fixed + n_random;
	//
	// mixed_object (no information)
	bool quasi_fixed   = false;
	bool bool_sparsity = true;
	CppAD::mixed::d_sparse_rcv A_rcv; // empty matrix
	cppad_mixed mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
	);
	//
	// jac_a1fun = f_u (theta, u)
	vector<a2_double> a2_theta_u(n_both);
	for(size_t j = 0; j < n_both; j++)
		a2_theta_u[j] = a2_double(0.0);
	CppAD::Independent(a2_theta_u);
	vector<a2_double> a2_jac(n_random);
	for(size_t j = 0; j < n_random; ++j)
	{	// Jacobian
		a2_double theta_j = a2_theta_u[j];
		a2_double u_j     = a2_theta_u[j + n_fixed];
		a2_double term    = ( u_j - theta_j * theta_j );
		a2_jac[j]         = a2_double(2) * term;
	}
	CppAD::ADFun<a1_double> jac_a1fun(a2_theta_u, a2_jac);
	//
	// ran_hes_uu_rc = sparsity pattern for f_{uu} (theta , u)
	CppAD::mixed::sparse_rc ran_hes_uu_rc(n_random, n_random, n_random);
	for(size_t k = 0; k < n_random; k++)
		ran_hes_uu_rc.set(k, k, k);
	//
	// ran_hes_uu_rcv = hessian f_{uu} (theta, u);
	CppAD::mixed::a1_sparse_rcv ran_hes_uu_rcv( ran_hes_uu_rc );
	for(size_t k = 0; k < n_random; ++k)
		ran_hes_uu_rcv.set(k, a1_double(2.0) );
	//
	// a1_ldlt_ran_hes
	CppAD::mixed::ldlt_eigen<a1_double> a1_ldlt_ran_hes(n_random);
	a1_ldlt_ran_hes.init( ran_hes_uu_rcv.pat() );
	a1_ldlt_ran_hes.update( ran_hes_uu_rcv );
	//
	// record W(beta, theta, u)
	vector<a1_double> a1_beta_theta_u(n_all);
	for(size_t j = 0; j < n_all; j++)
		a1_beta_theta_u[j] = 1.0;
	CppAD::Independent(a1_beta_theta_u);
	vector<a1_double> a1_u_opt = CppAD::mixed::order2random(
		mixed_object,
		n_fixed,
		n_random,
		jac_a1fun,
		a1_ldlt_ran_hes,
		a1_beta_theta_u
	);
	CppAD::ADFun<double> W(a1_beta_theta_u, a1_u_opt);
	//
	// beta_theta_u
	vector<double> beta_theta_u(n_all);
	for(size_t j = 0; j < n_fixed; ++j)
	{	// theta
		beta_theta_u[j + n_fixed]     = std::sqrt( double(j+1) );
		//
		// correspoding beta
		beta_theta_u[j]               = std::sqrt( double(j+1) );
		//
		// corresponding optimal u
		beta_theta_u[j + 2 * n_fixed] = double(j+1);
	}
	//
	// for each component of W
	for(size_t k = 0; k < n_random; k++)
	{	// compute the Hessian of k-th component of W w.r.t beta_theta_u
		vector<double> H = W.Hessian(beta_theta_u, k);
		//
		// Check Hessian w.r.t. beta is equal Hessian of optimal u w.r.t theta
		for(size_t i = 0; i < n_fixed; i++)
		{	for(size_t j = 0; j < n_fixed; j++)
			{	double H_ij = H[ i * n_all + j ];
				double check = 0.0;
				if( i == j && i == k )
					check = 2.0;
				ok &= CppAD::NearEqual(H_ij, check, eps, eps);
			}
		}
	}
	return ok;
}
// END C++
