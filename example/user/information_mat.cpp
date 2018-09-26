// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin information_mat.cpp$$
$spell
	CppAD
	cppad
	hes
	eval
	interp
	xam
$$

$section Observed Information Matrix: Example and Test$$

$head Model$$
$latex \[
	\B{p}( y_i | \theta , u ) \sim \B{N} ( u_i + \theta_0 , \theta_1^2 )
\] $$
$latex \[
	\B{p}( u_i | \theta ) \sim \B{N} ( 0 , 1 )
\] $$
$latex \[
	\B{p}( \theta ) \sim \B{N} ( 4 , 1 )
\] $$
It follows that the Laplace approximation is exact and
$latex \[
	\B{p}( y_i | \theta ) \sim \B{N} \left( \theta_0 , 1 + \theta_1^2 \right)
\] $$
The corresponding total objective for the fixed effects is
$latex \[
L( \theta ) =  C + \frac{1}{2} \left[
	( \theta_0 - 4 )^2 + ( \theta_1 - 4 )^2 +
	N \log \left( 1 + \theta_1^2 \right) +
	\left( 1 + \theta_1^2 \right)^{-1} \sum_{i=0}^{N-1} ( y_i - \theta_0 )^2
\right]
\] $$

The partial of the objective w.r.t $latex \theta_0$$ is
$latex \[
\partial_{\theta(0)} L ( \theta )
=
( \theta_0 - 4 )
-
\left( 1 + \theta_1^2 \right)^{-1} \sum_{i=0}^{N-1} ( y_i - \theta_0 )
\] $$

The partial of the objective w.r.t $latex \theta_1$$ is
$latex \[
\partial_{\theta(1)} L ( \theta )
=
( \theta_1 - 4 )
+
N \left( 1 + \theta_1^2 \right)^{-1} \theta_1
-
\left( 1 + \theta_1^2 \right)^{-2} \theta_1
	\sum_{i=0}^{N-1} ( y_i - \theta_0 )^2
\] $$

The second partial w.r.t $latex \theta_0$$, $latex \theta_0$$ is
$latex \[
\partial_{\theta(0)} \partial_{\theta(0)} L ( \theta )
=
1 + N \left( 1 + \theta_1^2 \right)^{-1}
\] $$

The second partial w.r.t $latex \theta_1$$, $latex \theta_1$$ is
$latex \[
\partial_{\theta(1)} \partial_{\theta(1)} L ( \theta )
=
1
+
N \left( 1 + \theta_1^2 \right)^{-1}
-
2 N \left( 1 + \theta_1^2 \right)^{-2} \theta_1^2
-
\left( 1 + \theta_1^2 \right)^{-2} \sum_{i=0}^{N-1} ( y_i - \theta_0 )^2
+
4 \left( 1 + \theta_1^2 \right)^{-3} \theta_1^2
	\sum_{i=0}^{N-1} ( y_i - \theta_0 )^2
\] $$

The second partial w.r.t $latex \theta_0$$, $latex \theta_1$$ is
$latex \[
\partial_{\theta(0)} \partial_{\theta(1)} L ( \theta )
=
2 \left( 1 + \theta_1^2 \right)^{-2} \theta_1
	\sum_{i=0}^{N-1} ( y_i - \theta_0 )
\] $$

$code
$srcfile%example/user/information_mat.cpp
	%0%// BEGIN C++%// END C++%1%$$
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
	using CppAD::mixed::s_vector;
	using CppAD::mixed::d_vector;
	using CppAD::mixed::a1_vector;
	//
	class mixed_derived : public cppad_mixed {
	private:
		const size_t          n_fixed_;
		const size_t          n_random_;
		const d_vector&       y_;
	// ----------------------------------------------------------------------
	public:
		// constructor
		mixed_derived(
			size_t                 n_fixed        ,
			size_t                 n_random       ,
			bool                   quasi_fixed    ,
			bool                   bool_sparsity  ,
			const d_sparse_rcv&    A_rcv          ,
			const d_vector&        y              ) :
			cppad_mixed(
				n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
			)                     ,
			n_fixed_(n_fixed)     ,
			n_random_(n_random)   ,
			y_(y)
		{	assert( n_fixed == 2);
			assert( y_.size() == n_random_ );
		}
	// ----------------------------------------------------------------------
		// implementation of ran_likelihood
		template <typename Vector>
		Vector template_ran_likelihood(
			const Vector&         theta  ,
			const Vector&         u      )
		{	typedef typename Vector::value_type scalar;

			assert( theta.size() == n_fixed_ );
			assert( u.size() == y_.size() );
			Vector vec(1);

			// initialize part of log-density that is always smooth
			vec[0] = scalar(0.0);

			// sqrt_2pi = CppAD::sqrt(8.0 * CppAD::atan(1.0) );

			for(size_t i = 0; i < n_random_; i++)
			{	scalar mu     = u[i] + theta[0];
				scalar sigma  = theta[1];
				scalar res    = (y_[i] - mu) / sigma;

				// p(y_i | u, theta)
				vec[0] += log(sigma) + res * res / scalar(2.0);
				// following term does not depend on fixed or random effects
				// vec[0] += log(sqrt_2pi);

				// p(u_i | theta)
				vec[0] += u[i] * u[i] / scalar(2.0);
				// following term does not depend on fixed or random effects
				// vec[0] += log(sqrt_2pi);
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
			const Vector&         fixed_vec  )
		{	typedef typename Vector::value_type scalar;

			assert( fixed_vec.size() == n_fixed_ );
			Vector vec(1);

			// initialize part of log-density that is smooth
			vec[0] = scalar(0.0);

			// compute these factors once
			scalar sqrt_2pi = scalar(
				 CppAD::sqrt( 8.0 * CppAD::atan(1.0)
			));

			for(size_t j = 0; j < n_fixed_; j++)
			{	scalar mu     = scalar(4.0);
				scalar sigma  = scalar(1.0);
				scalar res    = (fixed_vec[j] - mu) / sigma;

				// This is a Gaussian term, so entire density is smooth
				vec[0]  += res * res / scalar(2.0);
				// following term does not depend on fixed effects
				// vec[0]  += log(sqrt_2pi * sigma);
			}
			return vec;
		}
		// a1_vector version of fix_likelihood
		virtual a1_vector fix_likelihood(const a1_vector& fixed_vec)
		{	return template_fix_likelihood( fixed_vec ); }
	};
}

bool information_mat_xam(void)
{
	bool   ok = true;
	double inf = std::numeric_limits<double>::infinity();
	double eps = 10. * std::numeric_limits<double>::epsilon();
	//
	size_t n_data   = 10;
	size_t n_fixed  = 2;
	size_t n_random = n_data;
	d_vector
		fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
	fixed_lower[0] = - inf; fixed_in[0] = 2.0; fixed_upper[0] = inf;
	fixed_lower[1] = .01;   fixed_in[1] = 0.5; fixed_upper[1] = inf;
	//
	// explicit constriants (in addition to l1 terms)
	d_vector fix_constraint_lower(0), fix_constraint_upper(0);
	//
	d_vector data(n_data), random_in(n_random);
	for(size_t i = 0; i < n_data; i++)
	{	data[i]       = double(i + 1);
		random_in[i] = 0.0;
	}

	// object that is derived from cppad_mixed
	bool quasi_fixed   = false;
	bool bool_sparsity = false;
	d_sparse_rcv A_rcv; // empty matrix
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, data
	);
	mixed_object.initialize(fixed_in, random_in);

	// optimize the fixed effects using Newton method
	std::string fixed_ipopt_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           adaptive\n"
		"String  derivative_test_print_all yes\n"
		"Numeric tol                       1e-8\n"
		"Integer max_iter                  15\n"
	;
	std::string random_ipopt_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           second-order\n"
		"Numeric tol                       1e-8\n"
	;
	d_vector random_lower(n_random), random_upper(n_random);
	for(size_t i = 0; i < n_random; i++)
	{	random_lower[i] = -inf;
		random_upper[i] = +inf;
	}
	// optimize fixed effects
	d_vector fixed_scale = fixed_in;
	CppAD::mixed::fixed_solution solution = mixed_object.optimize_fixed(
		fixed_ipopt_options,
		random_ipopt_options,
		fixed_lower,
		fixed_upper,
		fix_constraint_lower,
		fix_constraint_upper,
		fixed_scale,
		fixed_in,
		random_lower,
		random_upper,
		random_in
	);
	d_vector random_opt = mixed_object.optimize_random(
		random_ipopt_options,
		solution.fixed_opt,
		random_lower,
		random_upper,
		random_in
	);
	// compute corresponding information matrix
	d_sparse_rcv
	information_rcv = mixed_object.information_mat(solution, random_opt);
	//
	// there are three non-zero entries in the lower triangle
	ok  &= information_rcv.nnz() == 3;
	//
	// theta
	const d_vector& theta( solution.fixed_opt );
	//
	// theta_1^2
	double theta_1_sq = theta[1] * theta[1];
	// 1 + theta_1^2
	double var     = 1.0 + theta_1_sq;
	// (1 + theta_1^2)^2
	double var_sq = var * var;
	// (1 + theta_1^2)^3
	double var_cube = var * var * var;
	//
	// sum of ( y_i - theta_0 ) and ( y_i - theta_0 )^2
	double sum    = 0.0;
	double sum_sq = 0.0;
	for(size_t i = 0; i < n_data; i++)
	{	sum    += data[i] - theta[0];
		sum_sq += (data[i] - theta[0]) * (data[i] - theta[0]);
	}
	//
	// check Hessian values
	const s_vector& row( information_rcv.row() );
	const s_vector& col( information_rcv.col() );
	const d_vector& val( information_rcv.val() );
	for(size_t k = 0; k < 3; k++)
	{	if( row[k] == 0 )
		{	// only returning lower triangle
			ok &= col[k] == 0;
			// value of second partial w.r.t. theta_0, theta_0
			double check = 1.0 + double(n_data) / var;
			ok &= CppAD::NearEqual(val[k], check, eps, eps);
		}
		else if( col[k] == 0 )
		{	// only other choice for row
			ok &= row[k] == 1;
			// value of second partial w.r.t. theta_0, theta_1
			double check = 2.0 * theta[1] * sum / var_sq;
			ok &= CppAD::NearEqual(val[k], check, eps, eps);
		}
		else
		{	// only case that is left
			ok &= row[k] == 1;
			ok &= col[k] == 1;
			// value of second partial w.r.t. theta_1, theta_1
			double check = 1.0 + double(n_data) / var;
			check       -= 2.0 * double(n_data) * theta_1_sq / var_sq;
			check       -= sum_sq / var_sq;
			check       += 4.0 * theta_1_sq * sum_sq / var_cube;
			ok &= CppAD::NearEqual(val[k], check, eps, eps);
		}
	}
	return ok;
}
// END C++
