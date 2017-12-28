// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-17 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin newton_step.cpp$$
$spell
	CppAD
	cppad
$$

$section newton_step: Example and Test$$

$head Private$$
This example is not part of the
$cref/cppad_mixed public API/public/$$.

$code
$srcfile%example/private/newton_step.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/newton_step.hpp>

bool newton_step_xam(void)
{
	bool   ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();
	using CppAD::vector;
	typedef CppAD::AD<double>    a1_double;
	typedef CppAD::AD<a1_double> a2_double;
	//
	// compute a1_ad_fun --------------------------------------------------
	size_t n_fixed  = 2;
	size_t n_random = 3;
	size_t n_both   = n_fixed + n_random;
	//
	// a1fun
	vector<a2_double> a2_theta_u(n_both);
	for(size_t j = 0; j < n_both; j++)
		a2_theta_u[j] = a2_double(0.0);
	CppAD::Independent(a2_theta_u);
	vector<a2_double> a2_f(1);
	a2_f[0] = a2_double(0.0);
	for(size_t j = 0; j < n_both; j++)
	{	// Hessian is diagonal matrix with 2 on diagonal
		a2_double term = a2_theta_u[j];
		a2_f[0]       += term * term;
	}
	CppAD::ADFun<a1_double> a1fun(a2_theta_u, a2_f);
	//
	// jac_a1fun
	CppAD::Independent(a2_theta_u);
	vector<a2_double> a2_jac(n_random);
	for(size_t j = 0; j < n_random; ++j)
	{	a2_double a2_uj = a2_theta_u[n_fixed + j];
		a2_jac[j]       = a2_double(2.0) * a2_uj;
	}
	CppAD::ADFun<a1_double> jac_a1fun(a2_theta_u, a2_jac);

	/// compute has_ran and a1_hes_ran
	//
	// Only requesting Hessian w.r.t. random effects so the pattern
	// w.r.t to fixed effects does not matter.
	CppAD::mixed::sparse_rc pattern(n_both, n_both, n_random);
	for(size_t k = 0; k < n_random; k++)
		pattern.set(k, n_fixed + k, n_fixed + k);
	//
	// structure that holds result of calculation
	CppAD::mixed::sparse_rcv    hes_ran( pattern );
	//
	// create newton_checkpoint ---------------------------------------------
	vector<double> theta(n_fixed), u(n_random);
	for(size_t j = 0; j < n_fixed; j++)
		theta[j] = 0.0;
	for(size_t j = 0; j < n_random; j++)
		u[j] = 0.0;
	//
	CppAD::mixed::newton_step newton_checkpoint;
	newton_checkpoint.initialize(
		a1fun, jac_a1fun, hes_ran, theta, u
	);
	//
	// use and test newton_checkpoint ---------------------------------------
	vector<a1_double> a1_theta_u_v(n_fixed + 2 * n_random);
	for(size_t j = 0; j < n_fixed + 2 * n_random; j++)
		a1_theta_u_v[j] = double(j);
	vector<a1_double> a1_logdet_step(1 + n_random);
	newton_checkpoint.eval(a1_theta_u_v, a1_logdet_step);
	//
	// check log of determinant
	double logdet = Value( a1_logdet_step[0] );
	double check  = log(2.0) * n_random;
	ok           &= CppAD::NearEqual(logdet, check, eps, eps);
	//
	// check the Newton step
	for(size_t j = 0; j < n_random; j++)
	{	double step  = Value( a1_logdet_step[1 + j] );
		double vj    = Value( a1_theta_u_v[n_fixed + n_random + j] );
		double check = vj / 2.0;
		ok          &= CppAD::NearEqual(step, check, eps, eps);
	}

	return ok;
}
// END C++
