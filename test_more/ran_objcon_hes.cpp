// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>

/*
f(theta, u)     = [ ( e^u theta - z )^2  + u^2 ] / 2

f_u             = (e^u theta - z) e^u theta + u

f_u,theta       = (e^u theta - z) e^u + (e^u)^2 theta

f_u,theta,theta = 2 (e^u)^2

f_u,u           = (e^u theta - z) e^u theta + (e^u theta)^2 + 1
                = 2 (e^u theta)^2 - z e^u theta + 1

f_u,u,theta     = 4 (e^u theta) e^u - z e^u

f_u,u,u         = 4 (e^u theta) e^u theta - z e^u theta


The Laplace approximation for the random part of the objective is
h(theta, u) = log[2(e^u theta)^2 - z e^u theta + 1] + f(theta, u) + constant

The optimal uhat(theta) satisfies:
	f_u[theta, uhat(theta)] = 0
Its derivative can be found by solving for uhat_theta(theta) in the following:
	0 = f_u,u[theta, uhat(theta)] uhat_theta(theta)
	  + f_u,theta[theta, uhat(theta)]
Its second derivative can be found by solving for uhat_theta,theta(theta) in
	0 = f_u,u[theta, uhat(theta)] uhat_theta,theta(theta)
	  + f_u,u,u[theta, uhat(theta)] [u_theta(theta)]^2
      + 2 f_u,u,theta[theta, uhat(theta)] uhat_theta(theta)
	  + f_u,theta,theta[theta, uhat(theta]
*/

namespace {
	using CppAD::vector;
	using CppAD::log;
	using CppAD::AD;
	using CppAD::mixed::sparse_rcv;
	//
	using CppAD::mixed::a1_double;
	using CppAD::mixed::a2_double;
	//
	template <class scalar>
	scalar neg_loglike(
		const scalar& theta ,
		const scalar& u     ,
		const scalar& z     )
	{	//
		// p(z|theta)
		scalar res = (z - exp(u) * theta );
		scalar sum = res * res / scalar(2.0);
		//
		// p(u|theta)
		sum += u * u / scalar(2.0);
		//
		return sum;
	}

	class mixed_derived : public cppad_mixed {
	private:
		const double z_;
	public:
		// constructor
		mixed_derived(
			size_t                               n_fixed       ,
			size_t                               n_random      ,
			bool                                 quasi_fixed   ,
			bool                                 bool_sparsity ,
			const CppAD::mixed::sparse_rcv&      A_rcv         ,
			double                               z             ) :
			cppad_mixed(
				n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
			),
			z_(z)
		{ }
		// implementation of ran_likelihood
		virtual vector<a2_double> ran_likelihood(
			const vector<a2_double>& theta  ,
			const vector<a2_double>& u      )
		{	vector<a2_double> result(1);
			result[0] =  neg_loglike(theta[0], u[0], a2_double(z_));
			return result;
		}
	public:
		//
	};
	// -----------------------------------------------------------------------
	// function used to check results
	void derivative(
		const double&  theta,
		const double&  u               ,
		const double&  z               ,
		double&        f_u             ,
		double&        f_u_theta       ,
		double&        f_u_theta_theta ,
		double&        f_u_u           ,
		double&        f_u_u_theta     ,
		double&        f_u_u_u         )
	{	double eu       = exp(u);
		double eutheta  = eu * theta;
		//
		f_u             = ( eutheta - z ) * eutheta + u;
		f_u_theta       = ( eutheta - z ) * eu + eu * eutheta;
		f_u_theta_theta = 2.0 * eu;
		f_u_u           = 2.0 * eutheta * eutheta - z * eutheta + 1;
		f_u_u_theta     = 4.0 * eutheta * eu - z * eu;
		f_u_u_u         = 4.0 * eutheta * eu * theta - z * eutheta;
		//
		return;
	};
};

bool ran_objcon_hes(void)
{	using std::fabs;
	bool   ok    = true;
	double inf   = std::numeric_limits<double>::infinity();
	double eps99 = std::numeric_limits<double>::epsilon() * 99.0;

	// object that is derived from cppad_mixed
	size_t n_fixed     = 1;
	size_t n_random    = 1;
	bool quasi_fixed   = false;
	bool bool_sparsity = true;
	CppAD::mixed::sparse_rcv A_rcv; // empty matrix
	double z           = 1.0;
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, z
	);

	// initialize
	vector<double> fixed_vec(1), random_vec(1);
	fixed_vec[0]  = 0.5;
	random_vec[0] = 0.0;
	mixed_object.initialize(fixed_vec, random_vec);

	// lower and upper limits for random effects
	vector<double> random_lower(n_random), random_upper(n_random);
	for(size_t i = 0; i < n_random; i++)
	{	random_lower[i] = -inf;
		random_upper[i] = +inf;
	}

	// optimize the random effects
	std::string options;
	options += "Integer print_level     0\n";
	options += "String  sb              yes\n";
	options += "String  derivative_test second-order\n";
	options += "Numeric tol             1e-11\n";
	vector<double> random_opt = mixed_object.optimize_random(
		options, fixed_vec, random_lower, random_upper, random_vec
	);

	// compute derivatives of the random likeilhood at (theta, uhat)
	double theta  = fixed_vec[0];
	double uhat   = random_opt[0];
	double f_u, f_u_theta, f_u_theta_theta, f_u_u, f_u_u_theta, f_u_u_u;
	derivative(
		// inputs
		theta           ,
		uhat            ,
		z               ,
		// oututs
		f_u             ,
		f_u_theta       ,
		f_u_theta_theta ,
		f_u_u           ,
		f_u_u_theta     ,
		f_u_u_u
	);
	// check optimality
	ok &= fabs( f_u ) < 1e-10;

	// record y = f(x) where x = (theta,u)
	vector<a1_double> ax(2), ay(1);
	ax[0] = theta;
	ax[1] = uhat;
	CppAD::Independent(ax);
	ay[0] = neg_loglike(ax[0], ax[1], a1_double(z));
	CppAD::ADFun<double> f(ax, ay);

	// check f_u
	vector<double> x1(2), y1(1);
	x1[0] = 0.0;
	x1[1] = 1.0;
	y1    = f.Forward(1, x1);
	ok   &= fabs( f_u - y1[0] ) < eps99;

	// check f_u_theta and f_u_u
	vector<double> w(1), dw(4);
	w[0] = 1.0;
	dw    = f.Reverse(2, w);
	ok   &= fabs( dw[0 * 2 + 1] - f_u_theta ) < eps99;
	ok   &= fabs( dw[1 * 2 + 1] - f_u_u ) < eps99;

	// recompute f_u_u
	vector<double> x2(2), y2(1);
	x2[0] = 0.0;
	x2[1] = 0.0;
	y2    = f.Forward(2, x2);
	ok   &= fabs( f_u_u / 2.0 - y2[0] ) < eps99;

	// check f_u_u_theta
	dw.resize(6);
	dw    = f.Reverse(3, w);
	ok   &= fabs( dw[ 0 * 3 + 2 ] - f_u_u_theta / 2.0 ) < eps99;
	ok   &= fabs( dw[ 1 * 3 + 2 ] - f_u_u_u / 2.0 ) < eps99;

	return ok;
}
