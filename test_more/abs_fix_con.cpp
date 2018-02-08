// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
// case with absolute values and fixed constraints

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
		double                sigma_;
		const vector<double>& z_;
	public:
		// constructor
		mixed_derived(
			size_t n_fixed                    ,
			size_t n_random                   ,
			bool   quasi_fixed                ,
			bool   bool_sparsity              ,
			const sparse_rcv&       A_rcv     ,
			double sigma                      ,
			const vector<double>& z           ) :
			cppad_mixed(
				n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
			),
			n_fixed_(n_fixed)                           ,
			sigma_(sigma)                               ,
			z_(z)
		{	assert(z.size() == n_fixed); }
		// implementation of fix_likelihood as p(z|theta) * p(theta)
		template <typename Vector>
		Vector template_fix_likelihood(
			const Vector& fixed_vec  )
		{	typedef typename Vector::value_type scalar;

			// initialize log-density
			Vector vec(1 + n_fixed_);
			vec[0] = scalar(0.0);

			// compute this factors once
			scalar sqrt_2 = scalar( CppAD::sqrt( 2.0 ) );

			for(size_t j = 0; j < n_fixed_; j++)
			{	// Data term
				scalar res   = z_[j] - CppAD::exp( fixed_vec[j] );
				res        /= scalar( sigma_ );
				vec[0]     += log(sqrt_2);
				vec[1 + j] += sqrt_2 * res;
			}
			return vec;
		}
		// a1_vector version of fix_likelihood
		virtual a1_vector fix_likelihood(const a1_vector& fixed_vec)
		{	return template_fix_likelihood( fixed_vec ); }
	public:
		// ------------------------------------------------------------------
		// User defined virtual functions
		//
		// fix_constraint
		template <typename Vector>
		Vector template_fix_constraint(
			const Vector& fixed_vec )
		{
			Vector ret_val(n_fixed_);
			ret_val = fixed_vec;
			return ret_val;
		}
		// a1_vector version of fix_constraint
		virtual a1_vector fix_constraint(const a1_vector& fixed_vec)
		{	return template_fix_constraint( fixed_vec ); }
		// ------------------------------------------------------------------
	};
}

bool abs_fix_con(void)
{
	bool   ok = true;
	double inf = std::numeric_limits<double>::infinity();
	double tol = 1e-8;
	//
	size_t n_fixed  = 3;
	vector<double> z(n_fixed);
	for(size_t i = 0; i < n_fixed; i++)
		z[i] = double(i+3);

	// fixed effects
	vector<double>
		fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed),
		fixed_constraint_lower(n_fixed), fixed_constraint_upper(n_fixed);
	for(size_t j = 0; j < n_fixed; j++)
	{	fixed_constraint_lower[j] = fixed_lower[j] = - inf;
		fixed_in[j]    = 0.0;
		fixed_constraint_upper[j] = fixed_upper[j] = inf;
	}
	//
	// set an upper limit on the first fixed effects
	// but do it implicitly though the constraint function
	fixed_constraint_upper[0] = CppAD::log( z[0] ) / 2.0 ;
	//
	// no random effects
	size_t n_random = 0;
	vector<double> random_in(0);
	vector<double> random_lower(n_random), random_upper(n_random);
	std::string random_ipopt_options = "";
	CppAD::mixed::sparse_rcv A_rcv; // empty matrix
	//
	// object that is derived from cppad_mixed
	bool quasi_fixed   = false;
	bool bool_sparsity = false;
	double sigma       = 1.0;
	mixed_derived mixed_object(
			n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, sigma, z
	);
	mixed_object.initialize(fixed_in, random_in);

	// optimize the fixed effects using quasi-Newton method
	std::string fixed_ipopt_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           adaptive\n"
		"Numeric tol                       1e-8\n"
		"Integer max_iter                  15\n"
	;
	vector<double> fixed_scale = fixed_in;
	CppAD::mixed::fixed_solution solution = mixed_object.optimize_fixed(
		fixed_ipopt_options,
		random_ipopt_options,
		fixed_lower,
		fixed_upper,
		fixed_constraint_lower,
		fixed_constraint_upper,
		fixed_scale,
		fixed_in,
		random_lower,
		random_upper,
		random_in
	);
	vector<double> fixed_out = solution.fixed_opt;
	//
	for(size_t j = 0; j < n_fixed; j++)
	{	double check = CppAD::log( z[j] );
		if( j == 0 )
			check = CppAD::log( z[j] ) / 2.0;
		ok &= fabs( fixed_out[j] - check ) <= tol;
	}

	return ok;
}
// END C++
