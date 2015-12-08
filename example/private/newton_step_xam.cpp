// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin newton_step_xam.cpp$$
$spell
	CppAD
	cppad
$$

$section C++ newton_step: Example and Test$$

$head Private$$
This example is not part of the
$cref/cppad_mixed public API/cppad_mixed_public/$$.

$code
$verbatim%example/private/newton_step_xam.cpp
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

	// create a function f(theta, u)
	size_t n_fixed  = 2;
	size_t n_random = 3;
	size_t n_both   = n_fixed + n_random;
	//
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
	CppAD::ADFun<a1_double> a1_adfun(a2_theta_u, a2_f);
	//
	vector<double> theta(n_fixed), u(n_random);
	for(size_t j = 0; j < n_fixed; j++)
		theta[j] = 0.0;
	for(size_t j = 0; j < n_random; j++)
		u[j] = 0.0;
	//
	CppAD::mixed::newton_step newton_atom;
	newton_atom.initialize(a1_adfun, theta, u);
	//
	vector<a1_double> a1_theta_u_v(n_fixed + 2 * n_random);
	for(size_t j = 0; j < n_fixed + 2 * n_random; j++)
		a1_theta_u_v[j] = double(j);
	vector<a1_double> a1_logdet_step(1 + n_random);
	newton_atom.eval(a1_theta_u_v, a1_logdet_step);
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
