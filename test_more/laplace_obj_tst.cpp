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
Test Laplace part of the the objective where non-linear and initial
random effects far from solution.
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

			// p(y | u, theta) , p(u_i | theta)
			scalar sigma  = theta[0];
			for(size_t i = 0; i < y_.size(); ++i)
			{	scalar res = (y_[i] - exp(u[i]) ) / sigma;
				vec[0] += log(sqrt_2pi * sigma) + res * res / scalar(2.0);
				vec[0] += log(sqrt_2pi) + u[i] * u[i] / scalar(2.0);
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

bool laplace_obj_tst(void)
{
	bool   ok = true;

	size_t n_data   = 2;
	size_t n_fixed  = 1;
	size_t n_random = n_data;
	d_vector y(n_data), theta(n_fixed), u(n_random);
	d_vector uhat(n_random);

	theta[0]  = 1.0;
	for(size_t i = 0; i < n_data; ++i)
	{	y[i] = double(10 * i);
		u[i] = 0.0;
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed   = false;
	bool bool_sparsity = true;
	CppAD::mixed::d_sparse_rcv A_rcv; // empty matrix
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, y
	);
	// Test not yet passing
	// mixed_object.initialize(theta, u);

	return ok;
}
// END C++
