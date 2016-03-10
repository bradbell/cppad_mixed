// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin capture_xam.cpp$$
$spell
	CppAD
	cppad
	Biometrics
	Royle
	Resample
	Poisson
	xam
	logit
	covariate
$$

$section A Capture Re-capture Example and Speed Test$$

$head Syntax$$
$codei%build/speed/capture_xam %random_seed%$$

$head random_seed$$
Is the random seed
$cref/s_in/manage_gsl_rng/new_gsl_rng/s_in/$$ used during the call to
$cref manage_gsl_rng$$.

$head Notation$$
$table
$latex T$$     $cnext number of sampling times
$rnext
$latex I$$     $cnext number of sampling locations
$rnext
$latex N_i$$   $cnext size of the population at $th i$$ location
$rnext
$latex y_{i,t}$$ $cnext
	number of captures at location $latex i$$ and time $latex t$$
$rnext
$latex M_i$$   $cnext
	maximum of the measured number of captures at $th i$$ location
$rnext
$latex p_{i,t}$$ $cnext
	probability of capture at location $latex i$$ and time $latex t$$
$rnext
$latex u_t$$   $cnext random effect for each sampling time
$rnext
$latex \theta_0$$ $cnext constant term in logit of probability
$rnext
$latex \theta_1$$ $cnext mean for $latex N_i$$ given $latex \theta$$
$rnext
$latex \theta_2$$ $cnext standard deviation of random effects $latex u_t$$
$tend

$head Count Data$$
We use a binomial distribution to model the
probability of $latex y_{i,t}$$ given $latex N_i$$ and $latex p_{i,t}$$; i.e,
$latex \[
\B{P} ( y_{i,t} | N_i , p_{i,t} )
=
\left( \begin{array}{c} N_i \\ y_{i,t} \end{array} \right)
p_{i,t}^{y(i,t)} \left( 1 - p_{i,t}^{y(i,t)} \right)
\] $$
Furthermore, we assume that this probability
is independent for each $latex (i, t)$$.

$head Capture Probability$$
Section 2.4 of the reference below suggest a
covariate model for the probability of capture.
We use a similar model defined by
$latex \[
	\R{logit} ( p_{i,t} ) = - u_t - \theta_0
\] $$
It follows that
$latex \[
	p_{i,t} = [ 1 + \exp(  u_t + \theta_0 ) ]^{-1}
\] $$

$head Population Probability$$
We use a Poisson distribution to model the
probability of $latex N_i$$ given $latex \theta_1$$; i.e.,
$latex \[
\B{P} ( N_i | \theta  )
=
\theta_1^{N(i)} \frac{ \exp[ - \theta_1 ] }{ N_i ! }
\] $$
We assume that this probability
is independent for each $latex i$$.

$head p(u|theta)$$
We use a normal distribution, with mean zero and standard deviation
$latex \theta_2$$,
for the distribution of the random effects $latex u$$
given the fixed effects $latex \theta$$; i.e.,
$latex \[
\B{p} ( u | \theta )
=
\prod_{t=0}^{T-1}
	\frac{1}{ \theta_2 \sqrt{ 2 \pi } }
		\exp \left[ - \frac{1}{2} \frac{ u_t^2 }{ \theta_2^2 } \right]
\] $$
This specifies the
$cref/random prior density
	/cppad_mixed
	/Notation
	/Random Prior Density, p(u|theta)
/$$.


$head p(y|theta,u)$$
We define the vector of maximum measurement for each location by
$latex \[
	M_i = \max \left\{ y_{i,0} , \cdots , y_{i, T-1} \right\}
\] $$
We define a probability of capture function for each location and time by
$latex \[
p_{i,t} ( \theta , u ) =  [ 1 + \exp( u_t + \theta_0 ) ]^{-1}
\] $$
The binomial likelihood for $latex y_{i, t}$$ given
$latex N_i$$, $latex \theta$$, and $latex u$$ is
$latex \[
B_{i,t} ( N_i, \theta , u )
=
\left( \begin{array}{c} {N(i)} \\ y_{i,t} \end{array} \right)
	p_{i,t} ( \theta , u)^{y(i,t)}
	\left( 1 - p_{i,t}( \theta , u)^{y(i,t)} \right)
\] $$
We do not know the population at each location $latex N_i$$,
but instead have a Poisson prior for $latex N_i$$.
The likelihood for the data at the $th i$$ location is
$latex \[
L_i ( \theta , u )
=
\sum_{k=M(i)}^{2 M(i) + 1}
	\theta_1^k \frac{ \exp( - \theta_1 ) }{ k ! }
		\prod_{t=0}^{T-1} B_{i,t} ( k , \theta , u )
\] $$
We assume that this is a good approximation; i.e.,
$latex \[
L_i ( \theta , u )
\approx
\sum_{k=M(i)}^{+\infty}
	\theta_1^k \frac{ \exp( - \theta_1 ) }{ k ! }
		\prod_{t=0}^{T-1} B_{i,t} ( k , \theta , u )
\] $$
Our model for the likelihood of the data at all the locations,
given the fixed and random effects, is
$latex \[
\B{p}( y | \theta , u )
=
\prod_{i=0}^{I-1} L_i ( \theta , u )
\] $$
This specifies the
$cref/random data density
	/cppad_mixed
	/Notation
	/Random Data Density, p(y|theta,u)
/$$.

$head A$$
The
$cref/random constraint
	/cppad_mixed
	/Problem
	/Random Constraints
/$$
for this example is
$latex \[
	0 = \hat{u}_0 ( \theta ) + \cdots + \hat{u}_{T-1} ( \theta )
\] $$
where $latex \hat{u} ( \theta )$$ is the
$cref/optimal random effects
	/cppad_mixed
	/Notation
	/Optimal Random Effects, u^(theta)
/$$.
The corresponding
$cref/random constraint matrix
	/cppad_mixed
	/Notation
	/Random Constraint Matrix, A
/$$
$latex A$$ is the row vector of size $latex T$$ with all ones; i.e.,
$latex \[
	A = [ 1 , \cdots , 1 ] \in \B{R}^{1 \times T}
\] $$

$head p(theta)$$
For this example there is no
$cref/fixed prior density
	/cppad_mixed
	/Notation
	/Fixed Prior Density, p(theta)
/$$
$latex \B{p}(\theta)$$.

$head p(z|theta)$$
For this example there is no
$cref/fixed data density
	/cppad_mixed
	/Notation
	/Fixed Data Density, p(z|theta)
/$$
$latex \B{p}(z | \theta)$$.

$head c(theta)$$
For this example there is no
$cref/fixed constraint function
	/cppad_mixed
	/Notation
	/Fixed Data Density, p(z|theta)
/$$
$latex \B{p}(z | \theta)$$.

$head Reference$$
J. Andrew Royle,
Biometrics 60, 108-115 March 2004,
$italic
N-Mixture Models for Estimating Population Size
from Spatially Replicated Counts.
$$

$code
$srcfile%speed/capture_xam.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$

$end
-----------------------------------------------------------------------------
*/
// BEGIN C++
# include <ctime>
# include <gsl/gsl_randist.h>
# include <cppad/utility.hpp> // CppAD::vector
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/manage_gsl_rng.hpp>
# include <cppad/mixed/configure.hpp>

namespace { // BEGIN_EMPTY_NAMESPACE

using CppAD::vector;
using std::exp;
using std::log;
using CppAD::mixed::sparse_mat_info;

// simulate data, y
void simulate(
	size_t                 I     ,
	size_t                 T     ,
	const vector<double>&  theta ,
	vector<size_t>&        y     )
{	assert( theta.size() == 3 );
	assert( y.size() == I * T );
	// random number generator
	gsl_rng* rng = CppAD::mixed::get_gsl_rng();
	//
	// simulate population sizes
	vector<double> N(I);
	double mu =  theta[1];
	for(size_t i = 0; i < I; i++)
		N[i] = gsl_ran_poisson(rng, mu );
	//
	// simulate random effects
	vector<double> u(T);
	double sigma = theta[2] * sqrt( double(T) / (double(T) - 1.0) );
	double sum   = 0.0;
	for(size_t t = 0; t < T; t++)
	{	u[t] = gsl_ran_gaussian(rng, sigma);
		sum += u[t];
	}
	// offset random effects to be mean zero (removes one degree of freedom)
	for(size_t t = 0; t < T; t++)
		u[t] = u[t] - sum / double(T);
	//
	// simulate data
	for(size_t i = 0; i < I; i++)
	{	for(size_t t = 0; t < T; t++)
		{	// probability of capture
			double ex = exp( u[t] + theta[0] );
			double p = 1.0 /( 1.0  + ex );
			y[ i * T + t ] = gsl_ran_binomial(rng, p, N[i]);
		}
	}
	//
	return;
}

// cppad_mixed derived class
class mixed_derived : public cppad_mixed {
private:
	const size_t          I_; // number of locations
	const size_t          T_; // number of times
	const vector<size_t>& y_; // reference to data values
	// -----------------------------------------------------------------
	// set by constructor and then effectively const
	vector<size_t>        M_;      // max number of captures at each location
	size_t                K_;      // practical bound on population size
	vector<double>        logfac_; // logfac_[k] = log( k! )
// ------------------------------------------------------------------------
public:
	// constructor
	mixed_derived(
		size_t                 I           ,
		size_t                 T           ,
		bool                   quasi_fixed ,
		const  sparse_mat_info& A_info     ,
		vector<size_t>&        y           )
		:
		// n_fixed = 3, n_random = T
		cppad_mixed(3, T, quasi_fixed, A_info) ,
		I_(I)            ,
		T_(T)            ,
		y_(y)
	{	// set M_ and K_
		M_.resize(I);
		K_ = 0;
		for(size_t i = 0; i < I; i++)
		{	M_[i] = 0;
			for(size_t t = 0; t < T; t++)
				M_[i] = std::max( M_[i], y[ i * T + t] );
			K_ = std::max( M_[i] , K_ );
		}
		K_ = 2 + 2 * K_ ;
		logfac_.resize(K_);
		logfac_[0]  = 0.0;
		logfac_[1]  = 0.0;
		for(size_t k = 2; k < K_; k++)
			logfac_[k] = log( double(k) ) + logfac_[k-1];
	}
	// implementaion of ran_likelihood
	template <class Float>
	vector<Float> implement_ran_likelihood(
		const vector<Float>&  theta  ,
		const vector<Float>&  u      )
	{	vector<Float> vec(1);
		Float one( 1.0 );
		Float two( 2.0 );
		Float pi2( 8.0 * std::atan(1.0) );
		Float eps = Float( 10.0 * std::numeric_limits<double>::epsilon() );
		//  ------------------------------------------------------------
		// log [ p(u | theta) ]
		//  ------------------------------------------------------------
		Float sig = theta[2];
		vec[0] -= Float(T_) * ( two * log(sig) + log( pi2 ) );
		for(size_t t = 0; t < T_; t++)
		{	Float w = u[t] / ( sig + eps );
			vec[0] -= w * w;
		}
		vec[0] /= two;
		//  ------------------------------------------------------------
		// log [ p(y | theta, u) ]
		//  ------------------------------------------------------------
		//
		// y_{i,t} * log( p_{i,t} ) and log( 1.0 - p_{i,t} )
		vector<Float> yit_log_p(I_ * T_), log_1p(I_ * T_);
		for(size_t i = 0; i < I_; i++)
		{	for(size_t t = 0; t < T_; t++)
			{	Float ex   = exp( u[t] + theta[0] );
				Float    p = one / (one + ex );
				yit_log_p[ i * T_ + t ]  = Float(y_[i*T_+t]) * log(p + eps);
				log_1p[i * T_ + t ]      = log(one - p + eps);
			}
		}
		//
		// log( theta[1]^k * exp( - theta[1] ) / k! )
		vector<Float> log_poisson(K_);
		for(size_t k = 0; k < K_; k++)
		{	log_poisson[k] = log( theta[1] + eps ) * Float(k);
			log_poisson[k] -= theta[1];
			log_poisson[k] -= logfac_[k];
		}
		//
		// loop over locations
		for(size_t i = 0; i < I_; i++)
		{	// initialize sum that defines L_i
			Float Li = Float(0.0);
			for(size_t k = M_[i]; k < 2 * M_[i] + 2; k++)
			{	// initialize log of term for this k
				//
				// terms that need to be calculated with Float
				Float float_sum = log_poisson[k];
				//
				// terms that dno not need to use Float
				double double_sum = 0.0;
				//
				// terms in likelihood for data given k that do not
				// depend on value of t
				double_sum     += T_ * logfac_[k];
				//
				// now compute terms that depend on t
				for(size_t t = 0; t < T_; t++)
				{	size_t yit = y_[ i * T_ + t ];
					//
					// log [ (k choose yit) / k! ]
					// where the k! term is added outside the loop
					double_sum += logfac_[yit] - logfac_[k - yit];
					// log [ pit^yit ]
					float_sum += yit_log_p[i * T_ + t];
					// log [ (1 - pit)^(k - yit) ]
					float_sum += Float(k - yit) * log_1p[i * T_ + t];
				}
				Li += exp( float_sum + double_sum );
			}
			vec[0] += log( Li );
		}
		vec[0] = - vec[0];
		return vec;
	}
// ------------------------------------------------------------------------
public:
	virtual vector<a2_double> ran_likelihood(
		const vector<a2_double>& fixed_vec  ,
		const vector<a2_double>& random_vec )
	{	return implement_ran_likelihood(fixed_vec, random_vec); }
	//
	virtual vector<a1_double> ran_likelihood(
		const vector<a1_double>& fixed_vec  ,
		const vector<a1_double>& random_vec )
	{	return implement_ran_likelihood(fixed_vec, random_vec); }
};
} // END_EMPTY_NAMESPACE

int main(int argc, char *argv[])
{	bool ok = true;
	using std::cout;
	using std::endl;
	//
	if( argc != 2 )
	{	std::cerr << "usage: " << argv[0] << " random_seed" << endl;
		std::exit(1);
	}
	//
	size_t random_seed = std::atoi( argv[1] );
	random_seed        = CppAD::mixed::new_gsl_rng( random_seed );
	//
	size_t n_fixed  = 3;
	std::time_t start_time = std::time( CPPAD_MIXED_NULL_PTR );
	// problem size
	size_t n_random = 30;
	size_t I = 25;
	//
	size_t T = n_random;
	vector<double> theta_sim(n_fixed);
	theta_sim[0] =   0.50;  // constant term in logit model
	theta_sim[1] =   5.0;   // mean population size
	theta_sim[2] =   1.00;  // standard deviation of random effects

	// simulate y
	vector<size_t> y(I * T);
	simulate(I, T, theta_sim, y);

	// lower and upper limits
	vector<double> fix_constraint_lower, fix_constraint_upper;
	vector<double>
		theta_lower(n_fixed), theta_in(n_fixed), theta_upper(n_fixed);
	// constant term
	theta_lower[0] = 0.0;
	theta_in[0]    =  theta_sim[0] / 2.0;
	theta_upper[0] = +2.0;
	// mean population
	theta_lower[1] = 1.0;
	theta_in[1]    = theta_sim[1] / 2.0;
	theta_upper[1] = 20.0;
	// standard deviation of random effects
	theta_lower[2] = 1e-2;
	theta_in[2]    = theta_sim[2] / 2.0;
	theta_upper[2] = 4.0;

	// constrain the sum of the random effects to be zero
	CppAD::mixed::sparse_mat_info A_info;
	A_info.resize(T);
	for(size_t t = 0; t < T; t++)
	{	A_info.row[t] = 0;
		A_info.col[t] = t;
		A_info.val[t] = 1.0;
	}

	// create derived object
	bool quasi_fixed = (random_seed % 2) == 0;
	mixed_derived mixed_object(I, T, quasi_fixed, A_info, y);

	// initialize point to start optimization at
	vector<double>  u_in(T);
	for(size_t t = 0; t < T; t++)
		u_in[t] = 0.0;
	std::map<std::string, size_t> size_map =
		mixed_object.initialize(theta_in, u_in);

	// print sizes
	cout << endl;
	std::map<std::string, size_t>::const_iterator itr;
	for(itr = size_map.begin(); itr != size_map.end(); ++itr)
		cout << itr->first << " = " << itr->second << endl;

	// optimize the fixed effects
	std::string fixed_options =
		"Integer print_level               5\n"
		"String  sb                        yes\n"
		"String  derivative_test           none\n"
		"String  derivative_test_print_all no\n"
		"Numeric tol                       1e-8\n"
		"Integer max_iter                  40\n"
	;
	std::string random_options =
		"Integer print_level 0\n"
		"String  sb          yes\n"
		"String  derivative_test none\n"
	;
	double inf = std::numeric_limits<double>::infinity();
	vector<double> u_lower(n_random), u_upper(n_random);
	for(size_t i = 0; i < n_random; i++)
	{	u_lower[i] = -inf;
		u_upper[i] = +inf;
	}
	// optimize fixed effects
	vector<double> theta_out = mixed_object.optimize_fixed(
		fixed_options,
		random_options,
		theta_lower,
		theta_upper,
		fix_constraint_lower,
		fix_constraint_upper,
		theta_in,
		u_lower,
		u_upper,
		u_in
	);
	// correspnding optimal random effects
	vector<double> u_out = mixed_object.optimize_random(
		random_options,
		theta_out,
		u_lower,
		u_upper,
		u_in
	);
	std::time_t end_time = std::time( CPPAD_MIXED_NULL_PTR );
	//
	// check random effects
	double sum = 0.0;
	for(size_t j = 0; j < n_random; j++)
		sum += u_out[j];
	ok &= std::fabs( sum ) < 1e-8;
	//
	// check reults
	for(size_t j = 0; j < n_fixed; j++)
		ok &= std::fabs( theta_out[j] / theta_sim[j] - 1.0 ) < 3e-1;
	//
	cout << "elapsed seconds = " << end_time - start_time << endl;
	cout << "random_seed = " << random_seed << endl;
	if( ok )
		cout << "capture_xam: OK" << endl;
	else
	{	cout << "capture_xam: Error" << endl;
		for(size_t j = 0; j < n_fixed; j++)
			cout << theta_out[j] / theta_sim[j] - 1.0 << endl;
	}
	//
	CppAD::mixed::free_gsl_rng();
	#
	if( ok )
		return 0;
	return 1;
}
// END C++
