// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */

# include <gsl/gsl_randist.h>
# include <cppad/utility.hpp> // CppAD::vector
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/manage_gsl_rng.hpp>

namespace { // BEGIN_EMPTY_NAMESPACE

using CppAD::vector;
using std::exp;
using std::log;
using std::endl;
using CppAD::mixed::sparse_rcv;

// simulate covariates, x, and data, y
void simulate(
	size_t                 N     ,
	size_t                 I     ,
	const vector<double>&  theta ,
	vector<size_t>&        y     )
{	assert( theta.size() == 1 );
	assert( y.size() == I );
	//
	gsl_rng* rng = CppAD::mixed::get_gsl_rng();
	//
	// simulate population sizes
	double p  = theta[0];
	for(size_t i = 0; i < I; i++)
		y[ i ] = gsl_ran_binomial(rng, p, N);
	return;
}

// cppad_mixed derived class
class mixed_derived : public cppad_mixed {
private:
	const size_t          N_;      // size of population
	const vector<size_t>& y_;      // reference to data values
	vector<double>        logfac_; // logfac_[k] = log( k ! )
// ------------------------------------------------------------------------
public:
	// constructor
	mixed_derived(
		size_t N                      ,
		vector<size_t>&  y            ,
		const CppAD::mixed::sparse_rcv&      A_rcv
	) :
		// n_fixed=1, n_random=0, quasi_fixed=false, bool_sparsity=true
		cppad_mixed(1, 0, false, true, A_rcv) ,
		N_(N)                                  ,
		y_(y)
	{	logfac_.resize(N+1);
		logfac_[0]  = 0.0;
		logfac_[1]  = 0.0;
		for(size_t k = 2; k <= N; k++)
			logfac_[k] = log( double(k) ) + logfac_[k-1];
	}
	// implementaion of ran_likelihood
	template <class Float>
	vector<Float> implement_fix_likelihood(const vector<Float>&  theta)
	{	vector<Float> vec(1);
		Float eps( 10.0 * std::numeric_limits<double>::epsilon() );
		//  ------------------------------------------------------------
		// log [ p(y | theta) ]
		//  ------------------------------------------------------------
		Float p      = theta[0];
		size_t I     = y_.size();
		Float sumlog = Float(0.0);
		Float term;
		for(size_t i = 0; i < I; i++)
		{	double yi = y_[i];
			// log (N choose yi)
			sumlog += Float( logfac_[N_] - logfac_[yi] - logfac_[N_ - yi] );
			// log ( p^yi )
			sumlog += Float( yi ) * log(p + eps);
			// log [ ( 1 - p^yi )^(N - yi) ]
			sumlog += Float(N_ - yi) * log( Float(1.0 + eps) - p );
		}
		vec[0] = - sumlog;
		return vec;
	}
// ------------------------------------------------------------------------
public:
	//
	virtual vector<a1_double> fix_likelihood(
		const vector<a1_double>& fixed_vec  )
	{	return implement_fix_likelihood<a1_double>(fixed_vec); }
};
} // END_EMPTY_NAMESPACE

bool binomial(void)
{	bool ok = true;
	size_t n_fixed = 1;
	size_t random_seed = CppAD::mixed::new_gsl_rng(0);

	// simulation parameters
	size_t N = 30;
	size_t I = 10 * 1000;
	vector<double> theta_sim(n_fixed);
	theta_sim[0] =   0.25;      // constant term in covariate model

	// simulate data
	vector<size_t> y(I);
	simulate(N, I, theta_sim, y);

	CppAD::mixed::sparse_rcv A_rcv; // empty matrix

	// create derived object
	mixed_derived mixed_object(N, y, A_rcv);

	// initialize point to start optimization at
	vector<double> theta_in( n_fixed ), u_in(0);
	for(size_t j = 0; j < n_fixed; j++)
		theta_in[j] = theta_sim[j];
	mixed_object.initialize(theta_in, u_in);

	// lower and upper limits
	vector<double> fix_constraint_lower, fix_constraint_upper;
	vector<double> theta_lower(n_fixed), theta_upper(n_fixed);
	// limits on p
	theta_lower[0] = 1e-4;
	theta_upper[0] = 1.0 - 1e-4;
	assert( theta_lower[0] <= theta_sim[0] );

	// optimize the fixed effects
	std::string fixed_ipopt_options =
		"Integer print_level               0\n"
		"String  sb                        yes\n"
		"String  derivative_test           second-order\n"
		"String  derivative_test_print_all yes\n"
		"Numeric tol                       1e-5\n"
		"Integer max_iter                  10\n"
	;
	std::string random_ipopt_options =
		"Integer print_level 0\n"
		"String  sb          yes\n"
		"String  derivative_test second-order\n"
	;
	vector<double> u_lower(0), u_upper(0);
	vector<double> theta_scale = theta_in;
	CppAD::mixed::fixed_solution solution = mixed_object.optimize_fixed(
		fixed_ipopt_options,
		random_ipopt_options,
		theta_lower,
		theta_upper,
		fix_constraint_lower,
		fix_constraint_upper,
		theta_scale,
		theta_in,
		u_lower,
		u_upper,
		u_in
	);
	vector<double> theta_out = solution.fixed_opt;
	//
	for(size_t j = 0; j < n_fixed; j++)
	{	// std::cout << endl << theta_out[j] / theta_sim[j] - 1.0 << endl;
		ok &= std::fabs( theta_out[j] / theta_sim[j] - 1.0) < 1e-2;
	}
	//
	if( ! ok )
		std::cout << endl << "random_seed = " << random_seed << endl;
	CppAD::mixed::free_gsl_rng();
	return ok;
}
