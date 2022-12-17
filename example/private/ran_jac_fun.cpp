/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin ran_jac_fun.cpp$$
$spell
	jac
	CppAD
	cppad
	interp
	xam
$$

$nospell
$bold This is old cppad_mixed documentation:$$ Here is a link to its
$href%http://bradbell.github.io/cppad_mixed%current documentation%$$.
$$
$section ran_jac_fun: Example and Test$$

$head Private$$
This example is not part of the
$cref/cppad_mixed public API/base_class/$$.

$code
$srcthisfile%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>

namespace {
	using CppAD::vector;
	using CppAD::log;
	using CppAD::AD;
	using CppAD::mixed::d_sparse_rcv;
	//
	using CppAD::mixed::a1_double;

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
			const CppAD::mixed::d_sparse_rcv&    A_rcv         ,
			const vector<double>&                y             ) :
			cppad_mixed(
				n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
			),
			y_(y)
		{ }
		// implementation of ran_likelihood
		template <typename Vector>
		Vector template_ran_likelihood(
			const Vector& theta  ,
			const Vector& u      )
		{	typedef typename Vector::value_type scalar;

			Vector vec(1);

			// initialize part of log-density that is always smooth
			vec[0] = scalar(0.0);

			for(size_t i = 0; i < y_.size(); i++)
			{	scalar mu     = u[i];
				scalar sigma  = theta[i];
				scalar res    = (y_[i] - mu) / sigma;

				// (do not need 2*pi inside of log)
				vec[0]  += (log(sigma) + res*res) / scalar(2.0);
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

bool ran_jac_fun_xam(void)
{
	bool   ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();
	using CppAD::mixed::a1_double;

	size_t n_data   = 10;
	size_t n_fixed  = n_data;
	size_t n_random = n_data;
	vector<double> data(n_data);
	vector<double> theta(n_fixed), u(n_random);

	for(size_t i = 0; i < n_data; i++)
	{	data[i]  = double((i + 1) * (i + 1) );
		theta[i] = std::sqrt( double(i + 1) );
		u[i]     = 0.0;
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed   = false;
	bool bool_sparsity = true;
	CppAD::mixed::d_sparse_rcv A_rcv; // empty matrix
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, data
	);
	mixed_object.initialize(theta, u);

	// compute Jacobian with respect to random effects
	vector<a1_double> a1_both(n_fixed + n_random);
	mixed_object.pack(theta, u, a1_both);
	vector<a1_double> a1_jac = mixed_object.ran_jac_a1fun_.Forward(0, a1_both);

	// check the Jacobian
	for(size_t i = 0; i < n_random; i++)
	{	double sigma  = theta[i];
		double mu     = u[i];
		double res    = (double(data[i]) - mu) / sigma;
		double check  = - res / sigma;
		//
		// std::cout << "jac[i] = " << jac[i];
		// std::cout << ", check = " << check << std::endl;
		ok              &= fabs( Value( a1_jac[i] ) / check - 1.0) <= eps;
	}

	return ok;
}
// END C++
