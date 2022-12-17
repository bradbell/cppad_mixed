/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin ran_like_hes.cpp$$
$spell
	CppAD
	cppad
	hes
$$

$nospell
$bold This is old cppad_mixed documentation:$$ Here is a link to its
$href%http://bradbell.github.io/cppad_mixed%current documentation%$$.
$$
$section ran_like_hes: Example and Test$$

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

bool ran_like_hes_xam(void)
{
	bool   ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();
	using CppAD::mixed::a1_double;
	using CppAD::mixed::a1_vector;

	size_t n_data   = 10;
	size_t n_fixed  = n_data;
	size_t n_random = n_data;
	vector<double> data(n_data);
	a1_vector      theta_u(n_fixed + n_random);
	vector<double> fixed_vec(n_fixed), random_vec(n_random);

	for(size_t i = 0; i < n_data; i++)
	{	data[i]   = double((i + 1) * (i + 1) );
		theta_u[i] = fixed_vec[i]  = std::sqrt( double(i + 1) );
		theta_u[i + n_fixed] = random_vec[i] = 0.0;
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed   = false;
	bool bool_sparsity = true;
	CppAD::mixed::d_sparse_rcv A_rcv; // empty matrix
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, data
	);
	mixed_object.initialize(fixed_vec, random_vec);

	// compute Jacobian with respect to random effects
	CppAD::mixed::a1_sparse_rcv ran_hes_uu_rcv = CppAD::mixed::ran_like_hes(
		n_fixed,
		n_random,
		mixed_object.ran_jac_a1fun_,
		mixed_object.a1_ldlt_ran_hes_.pattern(),
		theta_u
	);

	// check the Hessian
	CppAD::mixed::s_vector col_major = ran_hes_uu_rcv.col_major();
	size_t nnz = ran_hes_uu_rcv.nnz();
	ok &= nnz == n_random;
	for(size_t k = 0; k < nnz; k++)
	{	size_t i    = ran_hes_uu_rcv.row()[ col_major[k] ];
		size_t j    = ran_hes_uu_rcv.col()[ col_major[k] ];
		a1_double v = ran_hes_uu_rcv.val()[ col_major[k] ];
		ok &= i == k;
		ok &= j == k;
		a1_double sigma  = fixed_vec[i];
		a1_double check  = 1.0 / (sigma * sigma);
		//
		// std::cout << "hes[k] = " << v;
		// std::cout << ", check = " << check << std::endl;
		ok   &= fabs( v / check - 1.0) <= a1_double( eps );
	}
	//
	return ok;
}
// END C++
