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
	\B{p}( y_0 | \theta , u ) \sim \B{N} ( u_0 + \theta_0 , 1 )
\] $$
$latex \[
	\B{p}( u_0 | \theta ) \sim \B{N} ( 0 , 1 )
\] $$
It follows that the Laplace approximation is exact and
$latex \[
	\B{p}( y_0 | \theta ) \sim \B{N} ( \theta_0 , 2 )
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
	//
	class mixed_derived : public cppad_mixed {
	private:
		const vector<double>& y_;
	public:
		// constructor
		mixed_derived(
			size_t                               n_fixed       ,
			size_t                               n_random      ,
			bool                                 quasi_fixed   ,
			bool                                 bool_sparsity ,
			const CppAD::mixed::sparse_rcv&      A_rcv         ,
			const vector<double>&                y             ) :
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
			scalar sqrt_2pi = scalar(
				 CppAD::sqrt(8.0 * CppAD::atan(1.0)
			));

			for(size_t i = 0; i < y_.size(); i++)
			{	scalar mu     = u[i] + theta[0];
				scalar sigma  = scalar(1.0);
				scalar res    = (y_[i] - mu) / sigma;

				// p(y_i | u, theta)
				vec[0] += log(sqrt_2pi * sigma) + res*res / scalar(2.0);

				// p(u_i | theta)
				vec[0] += log(sqrt_2pi) + u[i] * u[i] / scalar(2.0);
			}
			return vec;
		}
		// a2_vector version of ran_likelihood
		// a3_vector version of ran_likelihood
		virtual a3_vector ran_likelihood(
			const a3_vector& fixed_vec, const a3_vector& random_vec
		)
		{	return template_ran_likelihood( fixed_vec, random_vec ); }
	};

	template <class Float>
	Float objective(const Float& y0, const Float& theta0)
	{
		Float sqrt_2pi = CppAD::sqrt(8.0 * CppAD::atan(1.0) );
		Float sigma    = 1.0;
		Float delta    = CppAD::sqrt( sigma * sigma + 1.0 );
		Float res      = (y0 - theta0) / delta;
		Float sum      = CppAD::log(sqrt_2pi * delta) + res*res / 2.0;

		return sum;
	}
}

bool laplace_obj_tst(void)
{
	bool   ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();

	size_t n_data   = 1;
	size_t n_fixed  = 1;
	size_t n_random = n_data;
	vector<double> y(n_data), theta(n_fixed), u(n_random);
	vector<double> uhat(n_random);

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

	// optimal random effects solve equation
	uhat[0] = (y[0] - theta[0]) / 2.0;

	// must factor f_{u,u} (theta, uhat)
	mixed_object.update_factor(theta, uhat);

	// compute random part of Laplace approximation using mixed_object
	double r_val = mixed_object.ran_obj_eval(theta, uhat);

	// tape laplace approximaiton for this case
	vector< AD<double> > a_theta(1), a_obj(1);
	a_theta[0] = theta[0];
	CppAD::Independent( a_theta );
	a_obj[0] = objective( AD<double>(y[0]), a_theta[0]);
	CppAD::ADFun<double> r(a_theta, a_obj);

	// check function value
	vector<double> check = r.Forward(0, theta);
	ok &= fabs( r_val / check[0] - 1.0 ) < eps;

	// check gradient value
	vector<double> r_fixed(n_fixed), d_theta(1);
	mixed_object.ran_obj_jac(theta, uhat, r_fixed);
	d_theta[0] = 1.0;
	check = r.Forward(1, d_theta);
	ok &= fabs( r_fixed[0] / check[0] - 1.0 ) < eps;

	// check hessian value
	vector<size_t> row_out, col_out;
	vector<double> val_out, weight(1);
	weight[0] = 1.0;
	mixed_object.laplace_obj_hes(
		theta, uhat, weight, row_out, col_out, val_out
	);
	check = r.Hessian(theta, 0);
	ok   &= val_out.size() == 1;
	ok   &= check.size() == 1;
	ok &= fabs( val_out[0] / check[0] - 1.0 ) < eps;

	return ok;
}
// END C++
