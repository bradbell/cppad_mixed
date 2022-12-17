/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin fix_like_hes.cpp$$
$spell
	CppAD
	cppad
	hes
	interp
	xam
$$

$nospell
$bold This is old cppad_mixed documentation:$$ Here is a link to its
$href%http://bradbell.github.io/cppad_mixed%current documentation%$$.
$$
$section fix_like_hes: Example and Test$$

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
		size_t                n_fixed_;
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
			n_fixed_(n_fixed) ,
			y_(y)
		{	assert( n_fixed == 2);
		}
		// implementation of ran_likelihood
		template <typename Vector>
		Vector template_ran_likelihood(
			const Vector& theta  ,
			const Vector& u      )
		{	typedef typename Vector::value_type scalar;

			Vector vec(1);

			// compute this factor once
			scalar sqrt_2pi = scalar(
				 CppAD::sqrt( 8.0 * CppAD::atan(1.0)
			));

			// initialize summation
			vec[0] = scalar(0.0);

			// for each data and random effect
			for(size_t i = 0; i < y_.size(); i++)
			{	scalar mu     = theta[0] + u[i];
				scalar sigma  = theta[1];
				scalar res    = (y_[i] - mu) / sigma;

				// This is a Gaussian term, so entire density is smooth
				vec[0]  += log(sqrt_2pi * sigma) + res * res / scalar(2.0);
			}
			return vec;
		}
		// a1_vector version of ran_likelihood
		virtual a1_vector ran_likelihood(
			const a1_vector& fixed_vec, const a1_vector& random_vec
		)
		{	return template_ran_likelihood( fixed_vec, random_vec ); }
		// implementation of fix_likelihood
		template <typename Vector>
		Vector template_fix_likelihood(
			const Vector& fixed_vec  )
		{	typedef typename Vector::value_type scalar;

			Vector vec(1);

			// initialize part of log-density that is smooth
			vec[0] = scalar(0.0);

			// compute these factors once
			scalar sqrt_2pi = scalar(
				 CppAD::sqrt( 8.0 * CppAD::atan(1.0)
			));

			for(size_t j = 0; j < n_fixed_; j++)
			{	scalar mu     = scalar(1.0);
				scalar sigma  = scalar(1.0);
				scalar res    = (fixed_vec[j] - mu);

				// This is a Gaussian term, so entire density is smooth
				vec[0]  += log(sqrt_2pi * sigma) + res * res / scalar(2.0);
			}
			return vec;
		}
		// a1_vector version of fix_likelihood
		virtual a1_vector fix_likelihood(const a1_vector& fixed_vec)
		{	return template_fix_likelihood( fixed_vec ); }
	public:
		//
	};
}

bool fix_like_hes_xam(void)
{
	bool   ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();

	size_t n_data   = 10;
	size_t n_fixed  = 2;
	size_t n_random = n_data;
	vector<double> data(n_data), fixed_vec(n_fixed), random_vec(n_random);

	fixed_vec[0] = 2.0;
	fixed_vec[1] = 0.5;
	for(size_t i = 0; i < n_data; i++)
	{	data[i]       = double(i + 1);
		random_vec[i] = double(i) / double(n_data);
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed   = false;
	bool bool_sparsity = false;
	CppAD::mixed::d_sparse_rcv A_rcv; // empty matrix
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, data
	);
	mixed_object.initialize(fixed_vec, random_vec);

	// compute prior hessian
	CppAD::vector<size_t> row, col;
	CppAD::vector<double> val, weight(1);
	weight[0] = 1.0;
	mixed_object.fix_like_hes(fixed_vec, weight, row, col, val);

	// initialize which variables have been found so far
	CppAD::vector<bool> found(3);
	for(size_t j = 0; j < n_fixed; j++)
		found[j] = false;

	// check derivatives
	for(size_t k = 0; k < row.size(); k++)
	{	size_t i = row[k];
		size_t j = col[k];
		double check = 1.0;
		ok      &= i == j;
		ok      &= ! found[i];
		ok      &= ( val[k] / check - 1.0) <= eps;
		found[i] = true;
	}
	for(size_t j = 0; j < n_fixed; j++)
		ok &= found[j];

	return ok;
}
// END C++
