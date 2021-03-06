// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-20 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin fix_likelihood.cpp$$
$spell
	CppAD
	cppad
	interp
	xam
$$

$section Random Likelihood: Example and Test$$

$code
$srcthisfile%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>


namespace {
	using CppAD::log;
	using CppAD::AD;
	//
	using CppAD::mixed::d_sparse_rcv;
	using CppAD::mixed::a1_double;
	using CppAD::mixed::d_vector;
	using CppAD::mixed::a1_vector;
	//
	class mixed_derived : public cppad_mixed {
	private:
		const d_vector&       z_;
	public:
		// constructor
		mixed_derived(
			size_t                 n_fixed       ,
			size_t                 n_random      ,
			bool                   quasi_fixed   ,
			bool                   bool_sparsity ,
			const d_sparse_rcv&    A_rcv         ,
			const d_vector&        z             ) :
			cppad_mixed(
				n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
			),
			z_(z)
		{ }
		// implementation of fix_likelihood
		template <typename Vector>
		Vector template_fix_likelihood(
			const Vector&         theta  )
		{	typedef typename Vector::value_type scalar;

			Vector vec(1);

			// compute this factor once
			double sqrt_2pi = CppAD::sqrt( 8.0 * CppAD::atan(1.0) );

			// initialize summation
			vec[0] = scalar(0.0);

			// for each data and random effect
			for(size_t i = 0; i < z_.size(); i++)
			{	scalar mu     = theta[i];
				scalar res    = (z_[i] - mu) / 1.0;

				// This is a Gaussian term, so entire density is smooth
				vec[0]  += log(sqrt_2pi) + res * res / scalar(2.0);
			}
			return vec;
		}
		// a1_vector version of fix_likelihood
		virtual a1_vector fix_likelihood(const a1_vector& fixed_vec)
		{	return template_fix_likelihood( fixed_vec ); }
	};
}

bool fix_likelihood_xam(void)
{
	bool   ok  = true;
	double pi  = 4.0 * std::atan(1.0);
	double eps = 100. * std::numeric_limits<double>::epsilon();
	//
	typedef cppad_mixed::a1_double a1_double;

	size_t n_data   = 10;
	size_t n_fixed  = n_data;
	size_t n_random = 0;
	d_vector    data(n_data);
	d_vector    fixed_vec(n_fixed), random_vec(n_random);
	a1_vector a1_fixed(n_fixed);

	for(size_t i = 0; i < n_data; i++)
	{	data[i]       = double(i + 1);
		//
		fixed_vec[i]  = 1.5;
		a1_fixed[i]   = a1_double( fixed_vec[i] );
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed   = true;
	bool bool_sparsity = true;
	d_sparse_rcv A_rcv; // empty matrix
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, data
	);
	mixed_object.initialize(fixed_vec, random_vec);

	// Evaluate fix_likelihood
	a1_vector a1_vec(1);
	a1_vec = mixed_object.fix_likelihood(a1_fixed);

	// check the random likelihood
	double sum = 0.0;
	for(size_t i = 0; i < n_data; i++)
	{	double mu     = fixed_vec[i];
		double res    = (data[i] - mu);
		sum          += (std::log(2 * pi) + res * res) / 2.0;
	}
	ok &= fabs( a1_vec[0] / a1_double(sum) - a1_double(1.0) ) < eps;

	return ok;
}
// END C++
