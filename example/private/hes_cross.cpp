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
$begin hes_cross.cpp$$
$spell
	CppAD
	cppad
	hes
	interp
	xam
$$

$section Hessian Cross Terms: Example and Test$$


$head Private$$
This example is not part of the
$cref/cppad_mixed public API/public/$$.

$code
$srcfile%example/private/hes_cross.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/pack.hpp>

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
			y_(y)
		{ }
		// implementation of ran_likelihood
		virtual vector<a2_double> ran_likelihood(
			const vector<a2_double>& theta  ,
			const vector<a2_double>& u      )
		{	vector<a2_double> vec(1);

			// initialize part of log-density that is always smooth
			vec[0] = a2_double(0.0);

			for(size_t i = 0; i < y_.size(); i++)
			{	a2_double mu     = u[i];
				a2_double sigma  = theta[i];
				a2_double res    = (y_[i] - mu) / sigma;

				// (do not need 2*pi inside of log)
				vec[0]  += (log(sigma) + res*res) / a2_double(2.0);
			}
			return vec;
		}
	public:
		//
		//
		virtual vector<a1_double> fix_likelihood(
			const vector<a1_double>& fixed_vec  )
		{	a1_vector vec(1);
			vec[0] = 0.0;
			return vec;
		}
	};
}

bool hes_cross_xam(void)
{
	bool   ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();

	size_t n_data   = 10;
	size_t n_fixed  = n_data;
	size_t n_random = n_data;
	vector<double> data(n_data);
	vector<double> theta(n_fixed), u(n_random);
	vector<double> fixed_vec(n_fixed), random_vec(n_random);

	for(size_t i = 0; i < n_data; i++)
	{	data[i]      = double( (i + 1) * (i + 1) );
		fixed_vec[i] = theta[i] = std::sqrt( double(i + 1) );
		random_vec[i] = u[i] = 0.0;
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed   = false;
	bool bool_sparsity = true;
	CppAD::mixed::sparse_rcv A_info; // empty matrix
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_info, data
	);
	mixed_object.initialize(theta, u);

	// number of non-zeros in Hessian cross terms
	ok &= mixed_object.hes_cross_.row.size() == n_random;
	ok &= mixed_object.hes_cross_.col.size() == n_random;

	// compute Hessian cross terms
	vector<double> w(1), val_out(n_random), both(n_fixed + n_random);
	w[0] = 1.0;
	CppAD::vector< std::set<size_t> > not_used;
	mixed_object.pack(fixed_vec, random_vec, both);
	mixed_object.ran_like_fun_.SparseHessian(
		both,
		w,
		not_used,
		mixed_object.hes_cross_.row,
		mixed_object.hes_cross_.col,
		val_out,
		mixed_object.hes_cross_.work
	);

	CppAD::vector<size_t>& row(mixed_object.hes_cross_.row);
	CppAD::vector<size_t>& col(mixed_object.hes_cross_.col);

	// check Hessian
	for(size_t k = 0; k < n_random; k++)
	{	ok     &= row[k] >= n_fixed;
		ok     &= col[k] < n_fixed;
		size_t i = row[k] - n_fixed;
		size_t j = col[k];
		ok      &= (i == j);
		//
		double sigma3 = fixed_vec[i] * fixed_vec[i] * fixed_vec[i];
		double check  = 2.0 * (data[i] - random_vec[i])  / sigma3;
		ok              &= fabs( val_out[k] / check - 1.0) <= eps;
	}
	for(size_t k = 1; k < n_random; k++)
	{	ok &= col[k-1] <= col[k];
		if( col[k-1] == col[k] )
			ok &= row[k-1] < row[k];
	}
	return ok;
}
// END C++
