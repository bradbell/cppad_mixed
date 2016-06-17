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
$begin ran_likelihood_hes_xam.cpp$$
$spell
	CppAD
	cppad
	interp
	xam
$$

$section Random Likelihood Hessian: Example and Test$$

$code
$srcfile%example/user/ran_likelihood_hes_xam.cpp
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
	using CppAD::mixed::sparse_mat_info;
	//
	typedef AD<double>    a1_double;
	typedef AD<a1_double> a2_double;

	class mixed_derived : public cppad_mixed {
	private:
		const vector<double>& y_;
	public:
		// constructor
		mixed_derived(
			size_t n_fixed                    ,
			size_t n_random                   ,
			bool   quasi_fixed                ,
			const  sparse_mat_info& A_info    ,
			const vector<double>& y           )
			:
			cppad_mixed(n_fixed, n_random, quasi_fixed, A_info) ,
			y_(y)
		{ }
		// ------------------------------------------------------------------
		// implementation of ran_likelihood
		virtual vector<a2_double> ran_likelihood(
			const vector<a2_double>& theta  ,
			const vector<a2_double>& u      )
		{	vector<a2_double> vec(1);

			// compute this factor once
			// sqrt_2pi = CppAD::sqrt( 8.0 * CppAD::atan(1.0) );

			// initialize summation
			vec[0] = a2_double(0.0);

			// for each data and random effect
			for(size_t i = 0; i < y_.size(); i++)
			{	a2_double mu     = exp( u[i] );
				a2_double sigma  = theta[i];
				a2_double res    = (y_[i] - mu) / sigma;

				// This is a Gaussian term, so entire density is smooth
				vec[0]  += log(sigma) + res * res / a2_double(2.0);
				// following term does not depend on fixed or random effects
				// vec[0]  += log(sqrt_2pi);
			}
			return vec;
		}
		// ------------------------------------------------------------------
		// ran_likelihood_hes
		vector<a1_double> ran_likelihood_hes(
			const vector<a1_double>& theta  ,
			const vector<a1_double>& u      ,
			const vector<size_t>&    row    ,
			const vector<size_t>&    col    )
		{	size_t K = row.size();
			assert( col.size() == K );

			// return value
			vector<a1_double> val(K);

			// for each component of the return value
			for(size_t k = 0; k < K; k++)
			{	// initialize it as zero
				val[k] = a1_double(0.0);

				// for this function, only the diagonal elements are non-zero
				if( row[k] == col[k] )
				{	size_t i = row[k];
					//
					a1_double mu        = exp( u[i] );
					a1_double sigma     = theta[i];
					a1_double res       = (y_[i] - mu) / sigma;
					a1_double res_ui    = - mu / sigma;
					a1_double res_ui_ui = - mu / sigma;
					a1_double sq_ui     = res * res_ui;
					a1_double sq_ui_ui  = res_ui * res_ui + res * res_ui_ui;
					val[k]              = sq_ui_ui;
				}
			}
			return val;
		}
		// ------------------------------------------------------------------
		// example a2 version of ran_likelihood
		// ------------------------------------------------------------------
		// example a1 version of ran_likelihood
	};
}

bool ran_likelihood_hes_xam(void)
{
	bool   ok  = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();
	typedef cppad_mixed::a1_double a1_double;

	size_t n_data   = 10;
	size_t n_fixed  = n_data;
	size_t n_random = n_data;
	vector<double>    data(n_data);
	vector<double>    fixed_vec(n_fixed), random_vec(n_random);
	vector<a1_double> a1_fixed(n_fixed), a1_random(n_random);
	vector<a2_double> a2_fixed(n_fixed), a2_random(n_random);

	for(size_t i = 0; i < n_data; i++)
	{	data[i]       = double(i + 1);
		//
		fixed_vec[i]  = 1.5;
		a1_fixed[i]   = fixed_vec[i];
		a2_fixed[i]   = a1_fixed[i];
		//
		random_vec[i] = 0.0;
		a1_random[i]  = random_vec[i];
		a2_random[i]  = a1_random[i];
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed = true;
	CppAD::mixed::sparse_mat_info A_info; // empty matrix
	mixed_derived mixed_object(n_fixed, n_random, quasi_fixed, A_info, data);
	mixed_object.initialize(fixed_vec, random_vec);

	// record Evaluation random likelihood
	CppAD::Independent(a2_random);
	vector<a2_double> a2_vec(1);
	a2_vec = mixed_object.ran_likelihood(a2_fixed, a2_random);
	CppAD::ADFun<a1_double> a1_f(a2_random, a2_vec);

	// sparsity pattern of Hessian for this function
	vector<size_t>    row(n_data), col(n_data);
	vector<a1_double> val(n_data);
	for(size_t i = 0; i < n_data; i++)
	{	row[i] = i;
		col[i] = i;
	};

	// Evaluate ran_likelihood_hes
	val = mixed_object.ran_likelihood_hes(a1_fixed, a1_random, row, col);

	// Compute the Hessian using AD
	vector<a1_double> w(1), check(n_random * n_random);
	w[0]  = 1.0;
	check = a1_f.Hessian(a1_random, w);

	// check the Hessian values
	// (when NDEBUG is not defined, cppad_mixed also does this check)
	for(size_t i = 0; i < n_data; i++)
	{	ok &= CppAD::NearEqual(
			Value(val[i]), Value(check[i * n_random + i]), eps, eps
		);
	}
	return ok;
}
// END C++
