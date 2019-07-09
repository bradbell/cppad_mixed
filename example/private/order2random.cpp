// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-19 University of Washington
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
$cref/cppad_mixed public API/base_class/$$.

$code
$srcfile%example/private/order2random.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/order2random.hpp>
# include <cppad/mixed/ran_like_hes.hpp>

namespace {
	using CppAD::vector;
	using CppAD::log;
	using CppAD::AD;
	using CppAD::mixed::d_sparse_rcv;
	//
	using CppAD::mixed::a1_double;
	using CppAD::mixed::a1_vector;

	class mixed_derived : public cppad_mixed {
	public:
		// constructor
		mixed_derived(
			size_t                               n_fixed       ,
			size_t                               n_random      ,
			bool                                 quasi_fixed   ,
			bool                                 bool_sparsity ,
			const CppAD::mixed::d_sparse_rcv&    A_rcv         ) :
		cppad_mixed(
			n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
		)
		{ }
		// implementation of ran_likelihood
		template <typename Vector>
		Vector template_ran_likelihood(
			const Vector& theta  ,
			const Vector& u      )
		{	typedef typename Vector::value_type scalar;
			assert( theta.size() == u.size() );

			Vector vec(1);

			// initialize part of log-density that is always smooth
			vec[0] = scalar(0.0);

			for(size_t j = 0; j < theta.size(); j++)
			{	scalar res   = exp(u[j]) - theta[j];
				vec[0]      += res * res;
			}
			return vec;
		}
		// a1_vector version of ran_likelihood
		virtual a1_vector ran_likelihood(
			const a1_vector& fixed_vec, const a1_vector& random_vec
		)
		{	return template_ran_likelihood( fixed_vec, random_vec ); }
	};
}

bool order2random_xam(void)
{
	bool   ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();
	//
	using CppAD::vector;
	typedef CppAD::AD<double>    a1_double;
	//
	// n_fixed
	size_t n_fixed     = 3;
	//
	// n_random
	size_t n_random    = n_fixed;
	//
	// n_both
	size_t n_both      = n_fixed + n_random;
	//
	// mixed_object (no information)
	bool quasi_fixed   = false;
	bool bool_sparsity = true;
	CppAD::mixed::d_sparse_rcv A_rcv; // empty matrix
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
	);
	vector<double> fixed_vec(n_fixed), random_vec(n_random);
	for(size_t j = 0; j < n_fixed; ++j)
	{	fixed_vec[j]  = 0.0;
		random_vec[j] = 0.0;
	}
	mixed_object.initialize(fixed_vec, random_vec);
	//
	// fun = f(theta, u)
	vector<a1_double> a1_theta_u(n_both), a1_theta(n_fixed), a1_u(n_random);
	for(size_t j = 0; j < n_both; j++)
		a1_theta_u[j] = 0.0;
	CppAD::Independent(a1_theta_u);
	for(size_t j = 0; j < n_fixed; ++j)
	{	a1_theta[j] = a1_theta_u[j];
		a1_u[j]     = a1_theta_u[j + n_fixed];
	}
	vector<a1_double> a1_f = mixed_object.ran_likelihood(a1_theta, a1_u);
	CppAD::ADFun<double> fun(a1_theta_u, a1_f);
	//
	// a1fun = f(theta, u)
	CppAD::ADFun<a1_double, double> a1fun;
	a1fun = fun.base2ad();
	//
	// jac_fun = f_u (theta, u)
	vector<a1_double> a1_fu(n_random);
	for(size_t j = 0; j < n_both; j++)
		a1_theta_u[j] = a1_double(0.0);
	CppAD::Independent(a1_theta_u);
	vector<a1_double> jac_all = a1fun.Jacobian(a1_theta_u);
	for(size_t j = 0; j < n_random; j++)
		a1_fu[j] = jac_all[j + n_fixed];
	CppAD::ADFun<double> jac_fun(a1_theta_u, a1_fu);
	//
	// jac_a1fun = f_u(theta, u)
	CppAD::ADFun<a1_double, double> jac_a1fun;
	jac_a1fun = jac_fun.base2ad();
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
	vector<a1_double> a1_beta(n_fixed);
	for(size_t j = 0; j < n_fixed; j++)
		a1_beta[j] = double(j + 1);
	for(size_t j = 0; j < n_both; j++)
		a1_theta_u[j] = double(j + 1);
	size_t abort_op_index = 0;
	bool   record_compare = true;
	CppAD::Independent(a1_beta, abort_op_index, record_compare, a1_theta_u);
	//
	// ran_hes_uu_rcv = f_uu(theta, u)
	ran_hes_uu_rcv = CppAD::mixed::ran_like_hes(
		n_fixed,
		n_random,
		jac_a1fun,
		a1_ldlt_ran_hes.pattern(),
		a1_theta_u
	);
	//
	// a1_ldlt_ran_hes_.update
	a1_ldlt_ran_hes.update( ran_hes_uu_rcv );
	//
	vector<a1_double> W = CppAD::mixed::order2random(
		n_fixed,
		n_random,
		jac_a1fun,
		a1_ldlt_ran_hes,
		a1_beta,
		a1_theta_u
	);
	CppAD::ADFun<double> Wfun(a1_beta, W);
	//
	// beta, theta_u
	vector<double> beta(n_fixed), theta_u(n_both);
	for(size_t j = 0; j < n_fixed; ++j)
	{	// theta[j]
		double theta_j = double(j+2);
		//
		// corresponding optimal u[j]
		double u_j     = std::log(theta_j);
		//
		// theta
		theta_u[j]  = theta_j;
		//
		// beta
		beta[j]     = theta_j;
		//
		// uhat
		theta_u[j + n_fixed] = u_j;
	}
	//
	// set dynamic parameters
	Wfun.new_dynamic(theta_u);
	//
	// for each component of W
	for(size_t k = 0; k < n_random; k++)
	{	// compute the Hessian of k-th component of W w.r.t beta
		vector<double> H = Wfun.Hessian(beta, k);
		//
		// Check Hessian w.r.t. beta is equal Hessian of optimal u w.r.t theta
		for(size_t i = 0; i < n_fixed; i++)
		{	for(size_t j = 0; j < n_fixed; j++)
			{	double H_ij = H[ i * n_fixed + j ];
				double check = 0.0;
				if( i == j && i == k )
				{	double theta_j = double(j+2);
					check =  - 1.0 / (theta_j * theta_j);
					// std::cout << "H_ij = " << H_ij << " , ";
					// std::cout << "check = " << check << "\n";
				}
				ok &= CppAD::NearEqual(H_ij, check, eps, eps);
			}
		}
	}
	return ok;
}
// END C++
