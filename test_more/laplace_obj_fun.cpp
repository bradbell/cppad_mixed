// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
Test Laplace part of the the objective: function, gradient, and Hessian
for a very simple case.

$latex \[
	\B{p}( y_0 | \theta , u ) \sim \B{N} ( u_0 , \theta_0^2 )
\] $$
$latex \[
	\B{p}( u_0 | \theta ) \sim \B{N} ( 0 , 1 )
\] $$
It follows that the Laplace approximation is exact and
$latex \[
	\B{p}( y_0 | \theta ) \sim \B{N} ( 0 , 1 + \theta_0^2 )
\] $$

*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/order2random.hpp>
# include <cppad/mixed/ran_like_hes.hpp>


namespace {
	using CppAD::vector;
	using CppAD::log;
	using CppAD::AD;
	using CppAD::mixed::d_sparse_rcv;
	using CppAD::mixed::d_vector;
	using CppAD::mixed::a1_double;
	using CppAD::mixed::a1_vector;
	using CppAD::mixed::sparse_rc;
	// -----------------------------------------------------------------------
	class mixed_derived : public cppad_mixed {
	private:
		const d_vector& y_;
	public:
		// constructor
		mixed_derived(
			size_t                               n_fixed       ,
			size_t                               n_random      ,
			bool                                 quasi_fixed   ,
			bool                                 bool_sparsity ,
			const CppAD::mixed::d_sparse_rcv&    A_rcv         ,
			const d_vector&                      y             ) :
			cppad_mixed(
				n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
			),
			y_(y)
		{	assert( n_fixed == 1);
		}
		// implementation of ran_likelihood
		template <typename Vector>
		Vector template_ran_likelihood(
			const Vector& theta  ,
			const Vector& u      )
		{	typedef typename Vector::value_type scalar;

			Vector vec(1);

			// initialize part of log-density that is always smooth
			vec[0] = scalar(0.0);

			// pi
			scalar sqrt_2pi = scalar( CppAD::sqrt(8.0 * CppAD::atan(1.0)) );

			// p(y_0 | u, theta)
			scalar sigma  = theta[0];
			scalar res    = (y_[0] - u[0]) / sigma;
			vec[0] += log(sqrt_2pi * sigma) + res * res / scalar(2.0);

			// p(u_i | theta)
			vec[0] += log(sqrt_2pi) + u[0] * u[0] / scalar(2.0);

			return vec;
		}
		// a1_vector version of ran_likelihood
		virtual a1_vector ran_likelihood(
			const a1_vector& fixed_vec, const a1_vector& random_vec
		)
		{	return template_ran_likelihood( fixed_vec, random_vec ); }
	};
	// -----------------------------------------------------------------------
	template <class Float>
	Float check_obj(const Float& y0, const Float& theta0)
	{
		Float sqrt_2pi = CppAD::sqrt(8.0 * CppAD::atan(1.0) );
		Float delta    = CppAD::sqrt( 1.0 + theta0 * theta0 );
		Float res      = y0 / delta;
		Float sum      = CppAD::log(sqrt_2pi * delta) + res * res / 2.0;

		return sum;
	}
	// -----------------------------------------------------------------------
	// special recording for debugging dynamic parameters

}

bool laplace_obj_fun(void)
{
	bool   ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();

	size_t n_data   = 1;
	size_t n_fixed  = 1;
	size_t n_random = n_data;
	d_vector y(n_data), theta(n_fixed), u(n_random);
	d_vector uhat(n_random);

	theta[0]  = 4.0;
	y[0]      = 1.0;
	u[0]      = 0.0;

	// object that is derived from cppad_mixed
	bool quasi_fixed   = false;
	bool bool_sparsity = true;
	CppAD::mixed::d_sparse_rcv A_rcv; // empty matrix
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, y
	);
	mixed_object.initialize(theta, u);

	// change value of theta
	theta[0] = 2.0;

	// optimal random effects solve equation
	// (uhat[0] - y[0]) / theta[0]^2 + uhat[0] = 0
	// uhat[0] ( 1 + theta[0]^2 ) = y[0]
	uhat[0] = y[0] / (1.0 + theta[0] * theta[0]);

	// must factor f_{u,u} (theta, uhat)
	mixed_object.update_factor(theta, uhat);

	// tape laplace approximaiton for this case
	vector<a1_double> a_theta(1), a_obj(1);
	a_theta[0] = theta[0];
	CppAD::Independent( a_theta );
	a_obj[0] = check_obj( a1_double(y[0]), a_theta[0]);
	CppAD::ADFun<double> r(a_theta, a_obj);

	// -----------------------------------------------------------------------
	// r(theta) using check_obj
	d_vector check = r.Forward(0, theta);
	//
	// r(theta) using ran_obj_eval
	double ran_obj = mixed_object.ran_obj_eval(theta, uhat);
	//
	// r(theta) using laplace_obj_fun
	d_vector theta_u(2);
	theta_u[0] = theta[0];
	theta_u[1] = uhat[0];
	mixed_object.laplace_obj_fun_.new_dynamic(theta_u);
	d_vector laplace_obj =
		mixed_object.laplace_obj_fun_.Forward(0, theta);
	//
	// check r(theta)
	ok &= fabs( ran_obj / check[0] - 1.0 ) < eps;
	ok &= fabs( laplace_obj[0] / check[0] - 1.0 ) < eps;
	//
	// -----------------------------------------------------------------------
	// Jacobian using check_obj
	check = r.Jacobian(theta);
	//
	// Jacobian using ran_obj_jac
	d_vector ran_jac(n_fixed), d_theta(1);
	mixed_object.ran_obj_jac(theta, uhat, ran_jac);
	//
	// Jacobian using laplace_obj_fun
	d_vector laplace_jac =
		mixed_object.laplace_obj_fun_.Jacobian(theta);
	//
	// check Jacobina
	ok &= fabs( ran_jac[0] / check[0] - 1.0 ) < eps;
	ok &= fabs( laplace_jac[0] / check[0] - 1.0 ) < eps;
	//
	// -----------------------------------------------------------------------
	// Hessian using check_obj
	check = r.Hessian(theta, 0);
	//
	// Hessian using laplace_obj_fun
	vector<size_t> row_out, col_out;
	d_vector val_out, weight(1);
	weight[0] = 1.0;
	mixed_object.laplace_obj_hes(
		theta, uhat, weight, row_out, col_out, val_out
	);
	ok   &= val_out.size() == 1;
	ok   &= check.size() == 1;
	ok &= fabs( val_out[0] / check[0] - 1.0 ) < eps;
	// =======================================================================
	// Test creating newton step version of objective with dynamic parameters
	// =======================================================================
	uhat[0] = y[0] / (1.0 + theta[0] * theta[0]);
	//
	// beta, theta_u, beta_theta_u
	a1_vector a1_beta(1), a1_theta_u(2), a1_beta_theta_u(3);
	a1_beta[0]         = theta[0];
	a1_theta_u[0]      = theta[0];
	a1_theta_u[1]      = uhat[0];
	a1_beta_theta_u[0] = theta[0];
	a1_beta_theta_u[1] = theta[0];
	a1_beta_theta_u[2] = uhat[0];
	//
	double pi   = CppAD::atan(1.0) * 4.0;
	double constant_term = CppAD::log(2.0 * pi) * double(n_random) / 2.0;
	// -----------------------------------------------------------------------
	// F: recording beta as independent, theta_u as dynamic parameters
	size_t abort_op_index = 0;
	bool   record_compare = true;
	CppAD::Independent(a1_beta, abort_op_index, record_compare, a1_theta_u);
	//
	// ran_hes_uu_rcv = f_uu(theta, u)
	CppAD::mixed::a1_sparse_rcv ran_hes_uu_rcv;
	ran_hes_uu_rcv = CppAD::mixed::ran_like_hes(
		n_fixed,
		n_random,
		mixed_object.ran_jac_a1fun_,
		mixed_object.a1_ldlt_ran_hes_.pattern(),
		a1_theta_u
	);
	//
	// a1_ldlt_ran_hes_.update
	mixed_object.a1_ldlt_ran_hes_.update( ran_hes_uu_rcv );
	//
	a1_vector a1_W = CppAD::mixed::order2random(
		n_fixed,
		n_random,
		mixed_object.ran_jac_a1fun_,
		mixed_object.a1_ldlt_ran_hes_,
		a1_beta,
		a1_theta_u
	);
	//
	// Evaluate random likelihood f(beta, W)
	a1_vector a1_beta_W(2);
	a1_beta_W[0] = a1_beta[0];
	a1_beta_W[1] = a1_W[0];
	a1_vector a1_F = mixed_object.ran_like_a1fun_.Forward(0, a1_beta_W);
	//
	// Evaluate log det f_{uu} ( beta , W )
	a1_vector a1_F_hes  = mixed_object.ran_like_a1fun_.Hessian(a1_beta_W, 0);
	ok                 &= a1_F_hes.size() == 4;
	a1_double a1_logdet = log( a1_F_hes[3] );
	//
	// create function object corresponding to (logdet f_uu) / 2 + f
	a1_F[0] += a1_logdet / 2.0 - constant_term;
	CppAD::ADFun<double> F(a1_beta, a1_F);
	// -----------------------------------------------------------------------
	// G: recording beta_theta_u as independent variables
	CppAD::Independent(a1_beta_theta_u);
	a1_beta[0]    = a1_beta_theta_u[0];
	a1_theta_u[0] = a1_beta_theta_u[1];
	a1_theta_u[1] = a1_beta_theta_u[2];
	//
	// ran_hes_uu_rcv = f_uu(theta, u)
	ran_hes_uu_rcv = CppAD::mixed::ran_like_hes(
		n_fixed,
		n_random,
		mixed_object.ran_jac_a1fun_,
		mixed_object.a1_ldlt_ran_hes_.pattern(),
		a1_theta_u
	);
	//
	// a1_ldlt_ran_hes_.update
	mixed_object.a1_ldlt_ran_hes_.update( ran_hes_uu_rcv );
	//
	a1_W = CppAD::mixed::order2random(
		n_fixed,
		n_random,
		mixed_object.ran_jac_a1fun_,
		mixed_object.a1_ldlt_ran_hes_,
		a1_beta,
		a1_theta_u
	);
	//
	// Evaluate random likelihood f(beta, W)
	a1_beta_W[0] = a1_beta_theta_u[0];
	a1_beta_W[1] = a1_W[0];
	a1_vector a1_G = mixed_object.ran_like_a1fun_.Forward(0, a1_beta_W);
	//
	// Evaluate log det f_{uu} (beta , W)
	a1_vector a1_G_hes = mixed_object.ran_like_a1fun_.Hessian(a1_beta_W, 0);
	ok                &= a1_G_hes.size() == 4;
	a1_logdet          = log( a1_G_hes[3] );
	//
	// create function object corresponding to (logdet f_uu) / 2 + f
	a1_G[0] += a1_logdet / 2.0 - constant_term;
	CppAD::ADFun<double> G(a1_beta_theta_u, a1_G);
	// -----------------------------------------------------------------------
	// Compare computations using F and G
	//
	// change value of the beta, theta and u
	d_vector beta(1), beta_theta_u(3);
	theta[0]        = theta[0] / 2.0;
	uhat[0]         = y[0] / (1.0 + theta[0] * theta[0]);
	beta[0]         = theta[0];
	theta_u[0]      = theta[0];
	theta_u[1]      = uhat[0];
	beta_theta_u[0] = theta[0];
	beta_theta_u[1] = theta[0];
	beta_theta_u[2] = uhat[0];
	F.new_dynamic(theta_u);
	//
	// function
	d_vector r_val = r.Forward(0, theta);
	d_vector F_val = F.Forward(0, beta);
	d_vector G_val = G.Forward(0, beta_theta_u);
	ok &= fabs( F_val[0] / r_val[0] - 1.0 ) < eps;
	ok &= fabs( G_val[0] / r_val[0] - 1.0 ) < eps;
	//
	// Jacobian
	d_vector r_jac = r.Jacobian(theta);
	d_vector F_jac = F.Jacobian(beta);
	d_vector G_jac = G.Jacobian(beta_theta_u);
	ok &= fabs( F_jac[0] / r_jac[0] - 1.0 ) < eps;
	ok &= fabs( G_jac[0] / r_jac[0] - 1.0 ) < eps;
	//
	// Hessian
	d_vector r_hes = r.Hessian(theta, 0);
	d_vector F_hes = F.Hessian(beta, 0);
	d_vector G_hes = G.Hessian(beta_theta_u, 0);
	// ok &= fabs( F_hes[0] / r_hes[0] - 1.0 ) < eps; // this check fails
	ok &= fabs( G_hes[0] / r_hes[0] - 1.0 ) < eps;
	ok &= fabs( F_hes[0] / r_hes[0] - 1.0 ) < eps;
	//
	// -----------------------------------------------------------------------
	// std::cout
	// << "val = " << r_val[0] << ", " << F_val[0] << ", " << G_val[0] << "\n"
	// << "jac = " << r_jac[0] << ", " << F_jac[0] << ", " << G_jac[0] << "\n"
	// << "hes = " << r_hes[0] << ", " << F_hes[0] << ", " << G_hes[0] << "\n";
	// -----------------------------------------------------------------------
	return ok;
}
// END C++
