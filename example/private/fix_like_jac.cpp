// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin fix_like_jac.cpp$$
$spell
	CppAD
	cppad
	jac
	interp
	xam
$$

$section fix_like_jac: Example and Test$$

$head Private$$
This example is not part of the
$cref/cppad_mixed public API/public/$$.

$code
$srcfile%example/private/fix_like_jac.cpp
	%0%// BEGIN C++%// END C++%1%$$
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
	using CppAD::mixed::sparse_rcv;
	//
	using CppAD::mixed::a1_double;
	using CppAD::mixed::a2_double;

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
			const CppAD::mixed::sparse_rcv&      A_info        ,
			const vector<double>&                y             ) :
			cppad_mixed(
				n_fixed, n_random, quasi_fixed, bool_sparsity, A_info
			),
			n_fixed_(n_fixed) ,
			y_(y)
		{	assert( n_fixed == 2);
		}
		// implementation of ran_likelihood
		virtual vector<a2_double> ran_likelihood(
			const vector<a2_double>& theta  ,
			const vector<a2_double>& u      )
		{	vector<a2_double> vec(1);

			// compute this factor once
			a2_double sqrt_2pi = a2_double(
				 CppAD::sqrt( 8.0 * CppAD::atan(1.0)
			));

			// initialize summation
			vec[0] = a2_double(0.0);

			// for each data and random effect
			for(size_t i = 0; i < y_.size(); i++)
			{	a2_double mu     = theta[0] + u[i];
				a2_double sigma  = theta[1];
				a2_double res    = (y_[i] - mu) / sigma;

				// This is a Gaussian term, so entire density is smooth
				vec[0]  += log(sqrt_2pi * sigma) + res * res / a2_double(2.0);
			}
			return vec;
		}
		// implementation of fix_likelihood
		virtual vector<a1_double> fix_likelihood(
			const vector<a1_double>& fixed_vec  )
		{	vector<a1_double> vec(1);

			// initialize part of log-density that is smooth
			vec[0] = a1_double(0.0);

			// compute these factors once
			a1_double mu     = a1_double(1.0);
			a1_double sqrt_2 = CppAD::sqrt( a1_double(2.0) );

			for(size_t j = 0; j < n_fixed_; j++)
			{
				// This is a Laplace term
				vec[0] += CppAD::log( sqrt_2 );

				// part of the density that needs absolute value
				vec.push_back(sqrt_2 * (fixed_vec[j] - mu) );
			}
			return vec;
		}
	public:
		//
	};
}

bool fix_like_jac_xam(void)
{
	bool   ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();
	double sqrt_2 = CppAD::sqrt(2.0);

	size_t n_data   = 10;
	size_t n_fixed  = 2;
	size_t n_random = n_data;
	vector<double> data(n_data), fixed_vec(n_fixed), random_vec(n_random);

	fixed_vec[0] = 2.0;
	fixed_vec[1] = 0.5;
	for(size_t i = 0; i < n_data; i++)
	{	data[i]       = double(i + 1);
		random_vec[i] = i / double(n_data);
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed   = false;
	bool bool_sparsity = false;
	CppAD::mixed::sparse_rcv A_info; // empty matrix
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_info, data
	);
	mixed_object.initialize(fixed_vec, random_vec);

	// compute prior jacobian
	CppAD::vector<size_t> row, col;
	CppAD::vector<double> val;
	mixed_object.fix_like_jac(fixed_vec, row, col, val);

	// initialize which rows have been found so far
	CppAD::vector<bool> found(3);
	for(size_t i = 0; i < 3; i++)
		found[i] = false;

	// check derivatives
	for(size_t k = 0; k < row.size(); k++)
	{	size_t i = row[k];
		size_t j = col[k];
		double check = sqrt_2;
		ok      &= i == j+1;
		ok      &= ! found[i];
		ok      &= ( val[k] / check - 1.0) <= eps;
		found[i] = true;
	}
	ok &= found[0] == false;
	ok &= found[1];
	ok &= found[2];

	return ok;
}
// END C++
