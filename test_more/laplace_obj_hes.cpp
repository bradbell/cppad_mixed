// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-17 University of Washington
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

f_u_u           = (e^u theta - z) e^u theta + (e^u theta)^2 + 1
                = 2 (e^u theta)^2 - z e^u theta + 1

The Laplace approximation for the random part of the objective is
h(theta, u) = log[ f_u_u(theta, uhat) ] / 2 + f(theta, u) + constant

The optimal uhat(theta) satisfies:
	f_u[theta, uhat(theta)] = 0

We define
u_0(theta) = uhat
u_1(theta) = u_0(theta)-f_u[theta, u_0(theta)]/f_u_u [theta, u_0(theta)]
u_2(theta) = u_1(theta)-f_u[theta, u_1(theta)]/f_u_u [theta, u_1(theta)]
*/

namespace {
	using CppAD::vector;
	using CppAD::log;
	using CppAD::AD;
	using CppAD::mixed::sparse_rcv;
	//
	using CppAD::mixed::a1_double;
	using CppAD::mixed::a2_double;
	using CppAD::mixed::a1_vector;
	using CppAD::mixed::a2_vector;
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
		template <typename Vector>
		Vector template_ran_likelihood(
			const Vector& theta  ,
			const Vector& u      )
		{	typedef typename Vector::value_type scalar;

			Vector result(1);
			result[0] =  neg_loglike(theta[0], u[0], scalar(z_));
			return result;
		}
		// a2_vector version of ran_likelihood
		// a3_vector version of ran_likelihood
		virtual a3_vector ran_likelihood(
			const a3_vector& fixed_vec, const a3_vector& random_vec
		)
		{	return template_ran_likelihood( fixed_vec, random_vec ); }
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
		double&        f_u_u           )
	{	double eu       = exp(u);
		double eutheta  = eu * theta;
		//
		f_u             = ( eutheta - z ) * eutheta + u;
		f_u_u           = 2.0 * eutheta * eutheta - z * eutheta + 1.0;
		//
		return;
	}
	// -----------------------------------------------------------------------
	template <class scalar>
	scalar u_2( const scalar& theta , const scalar& uhat  , const scalar& z)
	{
		scalar eu       = exp(uhat);
		scalar eutheta  = eu * theta;
		//
		scalar f_u       = ( eutheta - z ) * eutheta + uhat;
		scalar f_u_u     = 2.0 * eutheta * eutheta - z * eutheta + 1.0;
		scalar u_1       = uhat - f_u / f_u_u;
		//
		eu               = exp(u_1);
		eutheta          = eu * theta;
		f_u              = ( eutheta - z ) * eutheta + u_1;
		f_u_u            = 2.0 * eutheta * eutheta - z * eutheta + 1.0;
		scalar result    = u_1 - f_u / f_u_u;
		//
		return result;
	}
	// -----------------------------------------------------------------------
	template <class scalar>
	scalar h(const scalar& theta, const scalar& uhat, const scalar& z)
	{
		scalar f        = neg_loglike(theta, uhat, z);
		scalar eu       = exp(uhat);
		scalar eutheta  = eu * theta;
		scalar f_u_u    = 2.0 * eutheta * eutheta - z * eutheta + 1.0;
		scalar result   = log( f_u_u ) / 2.0 + f;
		return result;
	}
};

bool laplace_obj_hes(void)
{	using std::fabs;
	bool   ok    = true;
	double inf   = std::numeric_limits<double>::infinity();
	double eps99 = std::numeric_limits<double>::epsilon() * 99.0;

	// object that is derived from cppad_mixed
	size_t n_fixed     = 1;
	size_t n_random    = 1;
	bool quasi_fixed   = false;
	bool bool_sparsity = false;
	CppAD::mixed::sparse_rcv A_rcv; // empty matrix
	double z           = 2.0;
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, z
	);

	// initialize
	vector<double> fixed_vec(1), random_vec(1);
	fixed_vec[0]  = z;
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
	double f_u, f_u_u;
	derivative(
		// inputs
		theta           ,
		uhat            ,
		z               ,
		// outputs
		f_u             ,
		f_u_u
	);
	// check optimality
	ok &= fabs( f_u ) < 1e-10;

	// record y = f(x) where x = (theta,u)
	a1_vector ax(2), ay(1);
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

	// recompute f_u_u
	vector<double> x2(2), y2(1);
	x2[0] = 0.0;
	x2[1] = 0.0;
	y2    = f.Forward(2, x2);
	ok   &= fabs( f_u_u / 2.0 - y2[0] ) < eps99;

	// define reduced objective r(theta) = h[theta , u_2(theta) ]
	a1_vector atheta(1), ar(1);
	atheta[0] = theta;
	CppAD::Independent(atheta);
	a1_double au     = uhat;
	a1_double az     = z;
	a1_double au_2   = u_2(atheta[0], au, az);
	ar[0]            = h(atheta[0], au_2, az);
	CppAD::ADFun<double> r(atheta, ar);
	//
	// compute hessian both ways
	vector<double> weight(1), th(1), hes(1);
	th[0]     = theta;
	hes       = r.Hessian(th, 0);
	weight[0] = 1.0;
	vector<size_t> row, col;
	vector<double> val;
	mixed_object.laplace_obj_hes(fixed_vec, random_opt, weight, row, col, val);
	ok &= fabs( val[0] - hes[0] ) < eps99;
	//
	return ok;
}
