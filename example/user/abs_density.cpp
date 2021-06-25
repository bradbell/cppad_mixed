// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin abs_density.cpp$$
$spell
	CppAD
	cppad
	hes
	eval
	interp
	xam
$$

$section Absolute Value In Log-Density: Example and Test$$.

$head Model$$
$latex \[
	\B{p}( z_i | \theta ) \sim \B{L} ( \theta_i , \sigma )
\] $$
where $latex \B{L} ( \mu , \sigma )$$ is the Laplace distribution
with mean $latex \mu$$ and standard deviation $latex \sigma$$.
The corresponding fixed likelihood
$cref/g(theta)/theory/Fixed Likelihood, g(theta)/$$
is
$latex \[
g( \theta ) = \sum_{i} \left[
	\log ( \sigma \sqrt{2} )
	+
	\sqrt{2} \; \left| \frac{ z_i - \exp( \theta_i )}{\sigma} \right|
\right]
\] $$
The optimal solution, with no constraints and no prior on $latex \theta$$ is
$latex \[
	\hat{\theta}_i = \log( z_i )
\] $$

$code
$srcthisfile%0%// BEGIN C++%// END C++%1%$$
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
	using CppAD::mixed::d_vector;
	using CppAD::mixed::a1_vector;
	//
	class mixed_derived : public cppad_mixed {
	private:
		size_t                n_fixed_;
		double                sigma_;
		const d_vector&       z_;
	public:
		// constructor
		mixed_derived(
			size_t                 n_fixed       ,
			size_t                 n_random      ,
			bool                   quasi_fixed   ,
			bool                   bool_sparsity ,
			const d_sparse_rcv&    A_rcv         ,
			double                 sigma         ,
			const d_vector&        z             ) :
			cppad_mixed(
				n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
			),
			n_fixed_(n_fixed)  ,
			sigma_(sigma)      ,
			z_(z)
		{	assert(z.size() == n_fixed); }

		// template version of fix_likelihood; i.e., p(z|theta) * p(theta)
		template <typename Vector>
		Vector template_fix_likelihood(
			const Vector&         fixed_vec  )
		{	// scalar

			// initialize log-density
			Vector vec(1 + n_fixed_);
			vec[0] = 0.0;

			// compute this factors once
			a1_double sqrt_2 =  CppAD::sqrt( 2.0  );

			for(size_t j = 0; j < n_fixed_; j++)
			{	// Data term
				a1_double res   = z_[j] - CppAD::exp( fixed_vec[j] );
				res        /=  sigma_ ;
				// the following term does not depend on fixed effects
				// vec[0]     += log(sigma_ * sqrt_2);
				vec[1 + j] += sqrt_2 * res;
			}
			return vec;
		}
		// a1_vector version of fix_likelihood
		virtual a1_vector fix_likelihood(const a1_vector& fixed_vec)
		{	return template_fix_likelihood( fixed_vec ); }
	};
}

bool abs_density_xam(void)
{
	bool   ok = true;
	double inf = std::numeric_limits<double>::infinity();
	double tol = 1e-8;

	// fixed effects
	size_t n_fixed  = 3;
	d_vector
		fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
	for(size_t j = 0; j < n_fixed; j++)
	{	fixed_lower[j] = - inf;
		fixed_in[j]    = 0.0;
		fixed_upper[j] = inf;
	}
	//
	// no random effects
	size_t n_random = 0;
	d_vector random_in(0);
	d_vector random_lower(n_random), random_upper(n_random);
	std::string random_ipopt_options = "";
	//
	// no constriants
	d_vector fix_constraint_lower(0), fix_constraint_upper(0);
	//
	d_vector z(n_fixed);
	for(size_t i = 0; i < n_fixed; i++)
		z[i] = double(i+3);

	// object that is derived from cppad_mixed
	bool quasi_fixed   = false;
	bool bool_sparsity = true;
	d_sparse_rcv A_rcv; // empty matrix
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
		"String  derivative_test_print_all yes\n"
		"Numeric tol                       1e-8\n"
		"Integer max_iter                  15\n"
	;
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
	d_vector fixed_out = solution.fixed_opt;

	for(size_t j = 0; j < n_fixed; j++)
		ok &= fabs( fixed_out[j] - CppAD::log( z[j] ) ) <= tol;

	return ok;
}
// END C++
