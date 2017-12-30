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
$$

$section Random Likelihood Hessian: Example and Test$$

$code
$srcfile%example/user/ran_likelihood_hes.cpp
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
	using CppAD::mixed::sparse_rcv;
	using CppAD::mixed::d_vector;
	using CppAD::mixed::a1_vector;
	using CppAD::mixed::a3_vector;
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

			Vector vec(1);

			// initialize summation
			vec[0] = scalar(0.0);

			// for each data and random effect
			for(size_t i = 0; i < y_.size(); i++)
			{	scalar model  = theta[i]*exp( u[i] );
				scalar res    = y_[i] - model;

				// This is a Gaussian term, so entire density is smooth
				vec[0]  += res * res / scalar(2.0);

				// prior for random effects
				vec[0]  += u[i] * u[i] / scalar(2.0);
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

				// for this function, only the diagonal elements are non-zero
				assert( row[k] == col[k] );
				//
				size_t i = row[k];
				//
				scalar model     = theta[i]*exp( u[i] );
				scalar res       = y_[i] - model;
				//
				scalar res_ui    = - model;
				scalar res_ui_ui = - model;
				scalar sq_ui     = res * res_ui;
				scalar sq_ui_ui  = res_ui * res_ui + res * res_ui_ui;
				val[k]           = sq_ui_ui;
				//
				val[k]          += 1.0;
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
	size_t n_data   = 10;
	size_t n_fixed  = n_data;
	size_t n_random = n_data;
	d_vector data(n_data);
	d_vector fixed_lower(n_data), fixed_in(n_data), fixed_upper(n_data);
	d_vector random_lower(n_data), random_in(n_data), random_upper(n_data);
	for(size_t i = 0; i < n_data; i++)
	{	data[i]       = double(i + 1);
		//
		fixed_lower[i] = random_lower[i] = -inf;
		fixed_upper[i] = random_upper[i] = +inf;
		fixed_in[i]    = data[i];
		random_in[i]   = 0.0;
	}
	//
	// object that is derived from cppad_mixed
	bool quasi_fixed   = false;
	bool bool_sparsity = true;
	sparse_rcv A_rcv; // empty matrix
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, data
	);
	mixed_object.initialize(fixed_in, random_in);

	// optimize the fixed effects using quasi-Newton method
	std::string fixed_ipopt_options =
		"Integer print_level               5\n"
		"String  sb                        yes\n"
		"String  derivative_test           adaptive\n"
		"String  derivative_test_print_all yes\n"
		"Integer max_iter                  50\n"
		"Numeric tol                       1e-8\n"
	;
	std::string random_ipopt_options =
		"Integer print_level     5\n"
		"String  sb              yes\n"
		"String  derivative_test second-order\n"
		"Numeric tol             1e-10\n"
	;
	//
	for(size_t i = 0; i < n_random; i++)
	{	random_lower[i] = -inf;
		random_upper[i] = +inf;
	}
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

	std::cout << "fixed_out = " << fixed_out << "\n";
	for(size_t i = 0; i < n_data; i++)
		ok &= CppAD::NearEqual( data[i], fixed_out[i], 1e-7, 1e-7);
	ok &= mixed_object.ran_likelihood_hes_called();


	return ok;
}
// END C++
