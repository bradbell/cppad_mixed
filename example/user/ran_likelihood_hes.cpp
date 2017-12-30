// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-17 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin ran_likelihood_hes.cpp$$
$spell
	CppAD
	cppad
	interp
	xam
	exponentiating
$$

$section Random Likelihood Hessian: Example and Test$$

$head Exponentiating Random Effects$$
For this example there is,
on fixed effects $latex \theta_0$$ and  one data point $latex y_i$$ for each
random effect $latex u_i$$.
We use $latex e_i$$ to denote the measurement error for the $th i$$
measurement.
The simulated data is
$latex \[
	y_i = \theta_0 \exp [ u_i ] + e_i
\] $$
where $latex u_i \sim \B{N} ( 0 , \sigma_u )$$
and $latex e_i \sim \B{N} ( 0 , \sigma_y )$$.

$head Motivation$$
Note that for this model,
the Hessian with respect to the random effects is non-trivial and
easy to calculate.
Also not this is a good example of exponentiating
random effects in the data model.


$code
$srcfile%example/user/ran_likelihood_hes.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/manage_gsl_rng.hpp>
# include <gsl/gsl_randist.h>

namespace {
	using CppAD::log;
	using CppAD::AD;
	//
	using CppAD::mixed::sparse_rc;
	using CppAD::mixed::sparse_rcv;
	//
	using CppAD::mixed::d_vector;
	using CppAD::mixed::a1_vector;
	using CppAD::mixed::a3_vector;
	//
	double random_effects_std = 0.5;
	double measurement_std    = 0.7;
	//
	class mixed_derived : public cppad_mixed {
	private:
		const d_vector&       y_;
		bool                  ran_likelihood_hes_called_;
	public:
		// constructor
		mixed_derived(
			size_t                 n_fixed       ,
			size_t                 n_random      ,
			bool                   quasi_fixed   ,
			bool                   bool_sparsity ,
			const sparse_rcv&      A_rcv         ,
			const d_vector&        y             ) :
			cppad_mixed(
				n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
			),
			y_(y),
			ran_likelihood_hes_called_(false)
		{ }
		// ------------------------------------------------------------------
		// implementation of ran_likelihood
		template <typename Vector>
		Vector template_ran_likelihood(
			const Vector&         theta  ,
			const Vector&         u      )
		{	typedef typename Vector::value_type scalar;
			assert( theta.size() == 1 );

			Vector vec(1);

			// initialize summation
			vec[0] = scalar(0.0);

			// for each data and random effect
			for(size_t i = 0; i < y_.size(); i++)
			{	scalar model  = theta[0]*exp( u[i] );

				// a Gaussian
				scalar res = ( y_[i] - model ) / measurement_std;
				vec[0]    += res * res / scalar(2.0);

				// prior for random effects
				res      = u[i] / random_effects_std;
				vec[0]  += res * res / scalar(2.0);
			}
			return vec;
		}
		// a3_vector version of ran_likelihood
		virtual a3_vector ran_likelihood(
			const a3_vector& fixed_vec, const a3_vector& random_vec
		)
		{	return template_ran_likelihood( fixed_vec, random_vec ); }
		// ------------------------------------------------------------------
		// ran_likelihood_hes
		template <typename Vector>
		Vector template_ran_likelihood_hes(
			const Vector&            theta  ,
			const Vector&            u      ,
			const s_vector&          row    ,
			const s_vector&          col    )
		{	typedef typename Vector::value_type scalar;

			size_t K = row.size();
			assert( col.size() == K );

			// return value
			Vector val(K);

			// for each component of the return value
			for(size_t k = 0; k < K; k++)
			{	// initialize it as zero
				val[k] = scalar(0.0);
				assert( theta.size() == 1 );

				// for this function, only the diagonal elements are non-zero
				assert( row[k] == col[k] );
				//
				size_t i = row[k];
				//
				scalar model     = theta[0]*exp( u[i] );
				scalar res       = (y_[i] - model) / measurement_std;
				//
				scalar res_ui    = - model / measurement_std;
				scalar res_ui_ui = - model / measurement_std;
				scalar sq_ui     = res * res_ui;
				scalar sq_ui_ui  = res_ui * res_ui + res * res_ui_ui;
				val[k]           = sq_ui_ui;
				//
				double var       = random_effects_std * random_effects_std;
				val[k]          += 1.0 / var;
			}
			return val;
		}
		// a1_vector version of ran_likelihood_hes
		virtual a1_vector ran_likelihood_hes(
			const a1_vector& theta ,
			const a1_vector& u     ,
			const s_vector&  row   ,
			const s_vector&  col   )
		{	ran_likelihood_hes_called_ = true;
			return template_ran_likelihood_hes( theta, u, row, col );
		}
		bool ran_likelihood_hes_called(void) const
		{	return ran_likelihood_hes_called_; }
	};
}

bool ran_likelihood_hes_xam(void)
{
	bool   ok  = true;
	double inf = std::numeric_limits<double>::infinity();
	//
	// size_t random_seed = CppAD::mixed::new_gsl_rng(0);
	CppAD::mixed::new_gsl_rng(0);
	gsl_rng* rng = CppAD::mixed::get_gsl_rng();
	//
	size_t n_fixed   = 1;
	d_vector fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
	fixed_lower[0]   = 0.0;
	fixed_in[0]      = 5.0;
	fixed_upper[0]   = 10.0;
	//
	size_t n_data    = 1000;
	size_t n_random  = n_data;
	d_vector data(n_data);
	d_vector
		random_lower(n_random), random_in(n_random), random_upper(n_random);
	double true_fixed_0 = 2.0;
	for(size_t i = 0; i < n_data; i++)
	{	double u_i = gsl_ran_gaussian(rng, random_effects_std);
		double e_i = gsl_ran_gaussian(rng, measurement_std);
		data[i]    = true_fixed_0 * exp( u_i ) + e_i;
		//
		random_lower[i] = -inf;
		random_in[i]    = 0.0;
		random_upper[i] = +inf;
	}
	//
	bool quasi_fixed   = false;
	bool bool_sparsity = true;
	sparse_rcv A_rcv; // empty matrix
	mixed_derived mixed_object(
		n_fixed, n_data, quasi_fixed, bool_sparsity, A_rcv, data
	);
	mixed_object.initialize(fixed_in, random_in);

	// optimize the fixed effects using quasi-Newton method
	std::string fixed_ipopt_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           adaptive\n"
		"String  derivative_test_print_all yes\n"
		"Integer max_iter                  50\n"
		"Numeric tol                       1e-8\n"
	;
	std::string random_ipopt_options =
		"Integer print_level     0\n"
		"String  sb              yes\n"
		"String  derivative_test none\n"
		"Numeric tol             1e-10\n"
	;
	//
	// ------------------------------------------------------------------
	// optimize fixed effects
	d_vector fixed_scale = fixed_in;
	d_vector fix_constraint_lower(0), fix_constraint_upper(0);
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
	d_vector fixed_out = solution.fixed_opt;
	//
	ok &= CppAD::NearEqual(true_fixed_0, fixed_out[0], 1e-2, 1e-2);
	ok &= mixed_object.ran_likelihood_hes_called();
	//
	CppAD::mixed::free_gsl_rng();
	return ok;
}
// END C++
