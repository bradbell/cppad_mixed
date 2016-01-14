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
$begin fix_con_hes_xam.cpp$$
$spell
	CppAD
	cppad
	hes
	interp
	xam
$$

$section constraint_hes: Example and Test$$

$head Private$$
This example is not part of the
$cref/cppad_mixed public API/public/$$.

$code
$verbatim%example/private/fix_con_hes_xam.cpp
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

	class mixed_derived : public cppad_mixed {
	private:
		size_t                n_fixed_;
		const vector<double>& y_;
	public:
		// constructor
		mixed_derived(
			size_t n_fixed                    ,
			size_t n_random                   ,
			const vector<double>& y           ) :
			// quasi_fixed = false
			cppad_mixed(n_fixed, n_random, false) ,
			n_fixed_(n_fixed) ,
			y_(y)
		{	assert( n_fixed == 2);
		}
	private:
		// implementation of ran_likelihood
		template <class Float>
		vector<Float> implement_ran_likelihood(
			const vector<Float>& theta  ,
			const vector<Float>& u      )
		{	vector<Float> vec(1);

			// compute this factor once
			Float sqrt_2pi = Float( CppAD::sqrt( 8.0 * CppAD::atan(1.0) ) );

			// initialize summation
			vec[0] = Float(0.0);

			// for each data and random effect
			for(size_t i = 0; i < y_.size(); i++)
			{	Float mu     = theta[0] + u[i];
				Float sigma  = theta[1];
				Float res    = (y_[i] - mu) / sigma;

				// This is a Gaussian term, so entire density is smooth
				vec[0]  += log(sqrt_2pi * sigma) + res * res / Float(2.0);
			}
			return vec;
		}
		// implementation of fix_likelihood
		template <class Float>
		vector<Float> implement_fix_likelihood(
			const vector<Float>& fixed_vec  )
		{	vector<Float> vec(1);

			// initialize part of log-density that is smooth
			vec[0] = Float(0.0);

			// compute these factors once
			Float mu     = Float(1.0);
			Float sqrt_2 = CppAD::sqrt( Float(2.0) );

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
		virtual vector<a2_double> ran_likelihood(
			const vector<a2_double>& fixed_vec  ,
			const vector<a2_double>& random_vec )
		{	return implement_ran_likelihood(fixed_vec, random_vec); }
		virtual vector<a1_double> ran_likelihood(
			const vector<a1_double>& fixed_vec  ,
			const vector<a1_double>& random_vec )
		{	return implement_ran_likelihood(fixed_vec, random_vec); }
		//
		virtual vector<a1_double> fix_likelihood(
			const vector<a1_double>& fixed_vec  )
		{	return implement_fix_likelihood(fixed_vec); }
		//
		// constraint is 1/2 norm squared of the fixed effects
		virtual vector<a1_double> fix_constraint(
			const vector<a1_double>& fixed_vec  )
		{	assert( fixed_vec.size() == n_fixed_ );
			vector<a1_double> vec(1);
			vec[0] = 0.0;
			for(size_t j = 0; j < n_fixed_; j++)
				vec[0] += fixed_vec[j] * fixed_vec[j];
			vec[0] /= 2.0;
			return vec;
		}
	};
}

bool fix_con_hes_xam(void)
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
		random_vec[i] = i / double(n_data);
	}

	// object that is derived from cppad_mixed
	mixed_derived mixed_object(n_fixed, n_random, data);
	CppAD::mixed::sparse_mat_info A_info; // empty matrix
	mixed_object.initialize(A_info, fixed_vec, random_vec);

	// compute the constraint function and check result
	CppAD::vector<size_t> row, col;
	CppAD::vector<double> weight(1), val;
	weight[0] = 1.0;
	mixed_object.fix_con_hes(fixed_vec, weight, row, col, val);


	// check derivatives
	CppAD::vector<bool> found(n_fixed);
	for(size_t j = 0; j < n_fixed; j++)
		found[j] = false;
	ok &= row.size() == n_fixed;
	for(size_t k = 0; k < row.size(); k++)
	{	size_t i = row[k];
		size_t j = col[k];
		double check = 1.0;
		ok      &= i == j;
		ok      &= ( val[k] / check - 1.0) <= eps;
		ok      &= ! found[i];
		found[i] = true;
	}

	return ok;
}
// END C++
