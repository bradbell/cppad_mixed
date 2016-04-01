// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
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
	using CppAD::mixed::sparse_mat_info;

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
			const  sparse_mat_info& A_info    ,
			double sigma                      ,
			const vector<double>& z           ) :
			cppad_mixed(n_fixed, n_random, quasi_fixed, A_info) ,
			n_fixed_(n_fixed)                           ,
			sigma_(sigma)                               ,
			z_(z)
		{	assert(z.size() == n_fixed); }
	private:
		// implementation of fix_likelihood as p(z|theta) * p(theta)
		template <class Float>
		vector<Float> implement_fix_likelihood(
			const vector<Float>& fixed_vec  )
		{
			// initialize log-density
			vector<Float> vec(1 + n_fixed_);
			vec[0] = Float(0.0);

			// compute this factors once
			Float sqrt_2 = Float( CppAD::sqrt( 2.0 ) );

			for(size_t j = 0; j < n_fixed_; j++)
			{	// Data term
				Float res   = z_[j] - CppAD::exp( fixed_vec[j] );
				res        /= Float( sigma_ );
				vec[0]     += log(sqrt_2);
				vec[1 + j] += sqrt_2 * res;
			}
			return vec;
		}
	public:
		// ------------------------------------------------------------------
		// User defined virtual functions
		virtual vector<a1_double> fix_likelihood(
			const vector<a1_double>& fixed_vec  )
		{	return implement_fix_likelihood(fixed_vec); }
		//
		// fix_constraint
		virtual vector<a1_double> fix_constraint(
			const vector<a1_double>& fixed_vec )
		{	vector<a1_double> ret_val(n_fixed_);
			ret_val = fixed_vec;
			return ret_val;
		}
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
	std::string random_options = "";
	CppAD::mixed::sparse_mat_info A_info; // empty matrix
	//
	// object that is derived from cppad_mixed
	bool quasi_fixed = false;
	double sigma     = 1.0;
	mixed_derived mixed_object(
			n_fixed, n_random, quasi_fixed, A_info, sigma, z
	);
	mixed_object.initialize(fixed_in, random_in);

	// optimize the fixed effects using quasi-Newton method
	std::string fixed_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           second-order\n"
		"Numeric tol                       1e-8\n"
	;
	CppAD::mixed::fixed_solution solution = mixed_object.optimize_fixed(
		fixed_options,
		random_options,
		fixed_lower,
		fixed_upper,
		fixed_constraint_lower,
		fixed_constraint_upper,
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
		ok &= CppAD::abs( fixed_out[j] - check ) <= tol;
	}

	return ok;
}
// END C++
