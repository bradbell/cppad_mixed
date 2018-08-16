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

namespace {
	using CppAD::vector;
	using CppAD::log;
	using CppAD::AD;
	using CppAD::mixed::sparse_rcv;
	using CppAD::mixed::d_vector;
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
			const CppAD::mixed::sparse_rcv&      A_rcv         ,
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
		// a3_vector version of ran_likelihood
		virtual a3_vector ran_likelihood(
			const a3_vector& fixed_vec, const a3_vector& random_vec
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

bool laplace_obj_tst(void)
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
	CppAD::mixed::sparse_rcv A_rcv; // empty matrix
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
	vector< AD<double> > a_theta(1), a_obj(1);
	a_theta[0] = theta[0];
	CppAD::Independent( a_theta );
	a_obj[0] = check_obj( AD<double>(y[0]), a_theta[0]);
	CppAD::ADFun<double> r(a_theta, a_obj);

	// -----------------------------------------------------------------------
	// r(theta) using check_obj
	d_vector check = r.Forward(0, theta);
	//
	// r(theta) using ran_obj_eval
	double ran_obj = mixed_object.ran_obj_eval(theta, uhat);
	//
	// r(theta) using laplace_obj_fun
	d_vector beta_theta_u(3);
	beta_theta_u[0] = beta_theta_u[1] = theta[0];
	beta_theta_u[2] = uhat[0];
	d_vector laplace_obj =
		mixed_object.laplace_obj_fun_.Forward(0, beta_theta_u);
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
		mixed_object.laplace_obj_fun_.Jacobian(beta_theta_u);
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
	//
	return ok;
}
// END C++
