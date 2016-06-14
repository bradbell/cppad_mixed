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
$escape $$
$spell
	std
	CppAD
	cppad
	Biometrics
	Royle
	Resample
	Poisson
	xam
	logit
	covariate
	ipopt
	gsl_rng
$$

$section A Capture Re-capture Example and Speed Test$$

$head Syntax$$
$codei%build/speed/capture_xam  \
	%random_seed% \
	%number_fixed_samples% \
	%number_locations% \
	%number_times%  \
	%max_population% \
	%mean_population% \
	%mean_logit_probability% \
	%std_logit_probability% \
	%quasi_fixed% \
	%random_constraint% \
	%trace_optimize_fixed%
%$$

$head Reference$$
J. Andrew Royle,
Biometrics 60, 108-115 March 2004,
$italic
N-Mixture Models for Estimating Population Size
from Spatially Replicated Counts.
$$

$head random_seed$$
This is a non-negative integer equal to the
seed for the random number generator,
to be specific,
$cref/s_in/manage_gsl_rng/new_gsl_rng/s_in/$$ used during the call to
$code new_gsl_rng$$.

$subhead actual_seed$$
If $icode random_seed$$ is zero,
the system clock, instead of $icode random_seed$$,
is used to seed the random number generator.
The actual random seed $icode actual_seed$$ is printed
so that you can reproduce results when $icode random_seed$$ is zero.

$head number_fixed_samples$$
This is a positive integer equal to the number of samples simulated
from the posterior distribution for the fixed effects using
$cref sample_fixed$$.
The samples are used to approximation the
standard deviation for the optimal fixed effects.
One should use a large number of samples (at least 100)
for this purpose.

$subhead quasi_fixed$$
If $icode quasi_fixed$$ is true,
some initialization that was skipped during $cref initialize$$,
is done so that the Hessian of the total objective
$cref/L(theta)/theory/Objective/Total Objective, L(theta)/$$
can be computed.
In this case, the amount of memory used by the objects in the
$cref/mixed_derived/derived_ctor/mixed_derived/$$ object
to be the similar to when $icode quasi_fixed$$ is true.

$head number_locations$$
This is a positive integer equal to the
number of locations at which the measurements are made; i.e.
$latex R$$ in the reference.
Increasing this value increases the amount for each function evaluation,
but does not change the number of fixed or random effects.

$head number_times$$
This is a positive integer equal to the
number of times at which the measurements are made; i.e.,
$latex T$$ in the reference.
This is also equal to the number of random effects in the model.

$head max_population$$
This is a positive integer equal to the
maximum value in the finite summation with respect to
population size; i.e.,
$latex K$$ in the reference.
This must be greater than any of the simulated number of captures
at any location and time; i.e., and $latex y_{i,t}$$.
Also note that $icode max_population$$ does not affect the simulated
data $latex y_{i,t}$$.
A suggested value is five times $icode mean_population$$.
The value of $icode max_population$$ is large enough if increasing it
takes more time but does not make a difference in the optimal fixed effects
(use the same value for the other arguments and same actual seed).

$head mean_population$$
This is a positive floating point value equal to the
mean of the Poisson distribution for the population
used to simulate data values;
$latex \lambda$$ in reference.

$head mean_logit_probability$$
This is a positive floating point value equal to the
mean of the logit of the capture probability
(independent of the random effects)
used to simulate data values.
The is also equal to the mean of the random effects.

$head std_logit_probability$$
This is a positive floating point value equal to the
standard deviation of the logit of the capture probability
(independent of the random effects)
used to simulate data values.
The is also equal to the standard deviation of the random effects.

$head quasi_fixed$$
This is either $code true$$ or $code false$$.
If $icode quasi_fixed$$ is true,
it is also true in the $code cppad_mixed$$
$cref/derived class constructor/derived_ctor/quasi_fixed/$$.
In this case, the a quasi-Newton approximation for the Hessian
of the total objective
$cref/L(theta)/theory/Objective/Total Objective, L(theta)/$$.
Otherwise, the Hessian of the total objective is computed using the
approximate random objective
$cref/H(beta, theta, u)
	/theory
	/Approximate Random Objective, H(beta, theta, u)
/$$.

$head random_constraint$$
This is either $code true$$ or $code false$$.
If it is $code false$$, there is no
$cref/random constraint
	/cppad_mixed
	/Problem
	/Random Constraints
/$$
for this example.
If it is $code true$$,
the random constraint is
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

$head Example$$
The $code cppad_mixed$$ automated testing system uses the following
values for the arguments to $code capture_xam$$:
$code
$verbatim%speed/CMakeLists.txt
	%0%# BEGIN capture_xam arguments%# END capture_xam arguments%0%$$
$$

$head trace_optimize_fixed$$
This is either $code true$$ or $code false$$.
If it is true, a $icode%print_level% = 5%$$
$cref/trace/ipopt_trace/$$ of the fixed effects optimization
is included in the program output.
Otherwise the ipopt $icode print_level$$ is zero and
no such trace is printed.

$head Notation$$
$table
$latex R$$     $cnext
	number of sampling locations; i.e., $icode number_locations$$.
$rnext
$latex T$$     $cnext
	number of sampling times; i.e., $icode number_times$$.
$rnext
$latex K$$     $cnext
	maximum population in truncation of infinite summation; i.e.,
	$icode max_population$$
$rnext
$latex \theta_0$$
	$cnext mean for $latex N_i$$ given $latex \theta$$
	($latex \lambda$$ in reference); i.e.,
	$icode mean_population$$.
$rnext
$latex \theta_1$$
	$cnext mean of logit of capture probability; i.e.,
	$icode mean_logit_probability$$.
$rnext
$latex \theta_2$$
	$cnext standard deviation of logit of capture probability; i.e.,
	$icode std_logit_probability$$.
$rnext
$latex N_i$$   $cnext
	size of the population at $th i$$ location
$rnext
$latex y_{i,t}$$ $cnext
	number of captures at location $latex i$$ and time $latex t$$
	($latex n_{i,t}$$ in reference)
$rnext
$latex y_i$$ $cnext
	is the vector of captures at location $latex i$$
	$latex ( y_{i,0} , \ldots , y_{i, T-1} )$$.
$rnext
$latex M_i$$   $cnext
	maximum of captures at $th i$$ location
	($latex \max_t n_{i,t}$$ in reference)
$rnext
$latex q_t$$ $cnext
	capture probability at $latex t$$ same for all locations
	($latex p_{i,t}$$ in reference)
$rnext
$latex u_t$$
	$cnext random effect for each sampling time
$tend

$head p(y_it|N_i,q)$$
We use a binomial distribution to model the
probability of $latex y_{i,t}$$ given $latex N_i$$ and $latex q_t$$; i.e,
$latex \[
\B{p} ( y_{i,t} | N_i , q )
=
\left( \begin{array}{c} N_i \\ y_{i,t} \end{array} \right)
q_t^{y(i,t)} \left( 1 - q_t^{y(i,t)} \right)
\] $$
Furthermore, we assume that this probability
is independent for each $latex (i, t)$$.

$head q_t(theta,u)$$
Section 2.4 of the
$cref/reference/capture_xam.cpp/Reference/$$ suggests a
covariate model for the probability of capture.
We use a similar model defined by
$latex \[
	\R{logit} ( q_t ) = u_t + \theta_1
\] $$
It follows that
$latex \[
q_t( \theta , u) = [ 1 + \exp(- u_t - \theta_1 ) ]^{-1}
\] $$

$head p(N_i|theta)$$
We use a Poisson distribution to model the
probability of $latex N_i$$ given $latex \theta_0$$; i.e.,
$latex \[
\B{p} ( N_i | \theta  )
=
\theta_0^{N(i)} \frac{ \exp[ - \theta_0 ] }{ N_i ! }
\] $$
We assume these this probability
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
In $code cppad_mixed$$ notation, this specifies the
$cref/random prior density
	/cppad_mixed
	/Notation
	/Random Prior Density, p(u|theta)
/$$.

$head M_i$$
We define the vector of maximum measurement for each location by
$latex \[
	M_i = \max \left\{ y_{i,0} , \cdots , y_{i, T-1} \right\}
\] $$

$head p(y_i|theta,u)$$
The probability for $latex y_i$$ (the captures at location $latex i$$)
given $latex N_i$$, $latex \theta$$, and $latex u$$ is
$latex \[
\B{p}( y_i | N_i, \theta , u )
=
\prod_{t=0}^{T-1}
\left( \begin{array}{c} {N(i)} \\ y_{i,t} \end{array} \right)
	q_t ( \theta , u)^{y(i,t)}
	\left( 1 - q_t( \theta , u)^{y(i,t)} \right)
\] $$
We do not know the population at each location $latex N_i$$,
but instead have a Poisson prior for $latex N_i$$.
We sum with respect to the possible values for $latex N_i$$
to get the probability of $latex y_i$$ given
$latex \theta$$ and $latex u$$.
$latex \[
\B{p}( y_i | \theta , u )
=
\sum_{k=0}^K \B{p}( y_i | N_i=k, \theta , u ) \B{p}( N_i=k | \theta )
\] $$
where $latex k$$ is the possible values for $latex N_i$$.
Note that $latex K$$ should be plus infinity, but we use a fixed
finite value for $latex K$$ as an approximation for the infinite sum.

$head p(y|theta,u)$$
Our model for the likelihood of the data at all the locations,
given the fixed and random effects, is
$latex \[
\B{p}( y | \theta , u )
=
\prod_{i=0}^{R-1} \B{p}( y_i | \theta , u )
\] $$
Expressed in terms of fixed effects $latex \theta$$,
the random effects $latex u$$,
and the data $latex y$$, this is
$latex \[
\B{p}( y | \theta , u )
=
\prod_{i=0}^{R-1}
\left[
\sum_{k=0}^{K-1}
\theta_0^k \frac{ \exp[ - \theta_0 ] }{ k ! }
\prod_{t=0}^{T-1}
\left( \begin{array}{c} {k} \\ y_{i,t} \end{array} \right)
	q_t ( \theta , u)^{y(i,t)}
	\left( 1 - q_t( \theta , u)^{y(i,t)} \right)
\right]
\] $$
In $code cppad_mixed$$ notation, this specifies the
$cref/random data density
	/cppad_mixed
	/Notation
	/Random Data Density, p(y|theta,u)
/$$.

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
	bool                   random_constraint ,
	size_t                 R                 ,
	size_t                 T                 ,
	const vector<double>&  theta             ,
	vector<size_t>&        y                 )
{	assert( theta.size() == 3 );
	assert( y.size() == R * T );
	//
	// random number generator
	gsl_rng* rng = CppAD::mixed::get_gsl_rng();
	//
	// simulate population sizes
	vector<double> N(R);
	double mu =  theta[0];
	for(size_t i = 0; i < R; i++)
		N[i] = gsl_ran_poisson(rng, mu );
	//
	// simulate random effects
	vector<double> u(T);
	double sigma = theta[2];
	double sum   = 0.0;
	for(size_t t = 0; t < T; t++)
	{	u[t] = gsl_ran_gaussian(rng, sigma);
		sum += u[t];
	}
	// adjust the random effects when using random constraint
	if( random_constraint )
	{	for(size_t t = 0; t < T; t++)
		{	// simulate mean zero values
			u[t] = u[t] - sum / double(T);
			//
			// correct for loss of one degree of freedom
			u[t] = u[t] * sqrt( double(T) / double(T-1) );
		}
	}
	// simulate data
	for(size_t i = 0; i < R; i++)
	{	for(size_t t = 0; t < T; t++)
		{	// probability of capture
			double ex = exp( - u[t] - theta[1] );
			double q = 1.0 /( 1.0  + ex );
			y[ i * T + t ] = gsl_ran_binomial(rng, q, N[i]);
		}
	}
	//
	return;
}

// cppad_mixed derived class
class mixed_derived : public cppad_mixed {
private:
	const size_t          R_; // number of locations
	const size_t          T_; // number of times
	size_t                K_; // max used when summing over population size
	const vector<size_t>& y_; // reference to data values
	// -----------------------------------------------------------------
	// set by constructor and then effectively const
	vector<size_t>  M_;       // max number of captures at each location
	vector<double>  logfac_;  // logfac_[k] = log( k! )
	vector<double>  log_pik_; // used by implement_ran_likelihood
// ------------------------------------------------------------------------
public:
	// constructor
	mixed_derived(
		size_t                  R           ,
		size_t                  T           ,
		size_t                  K           ,
		bool                    quasi_fixed ,
		const  sparse_mat_info& A_info     ,
		vector<size_t>&         y           ,
		vector<double>&         fixed_in    ,
		vector<double>&         random_in   )
		:
		// n_fixed = 3, n_random = T
		cppad_mixed(3, T, quasi_fixed, A_info) ,
		R_(R)            ,
		T_(T)            ,
		K_(K)            ,
		y_(y)
	{	assert( fixed_in.size() == 3 );
		assert( random_in.size() == T );
		//
		// set M_
		M_.resize(R);
		for(size_t i = 0; i < R; i++)
		{	M_[i] = 0;
			for(size_t t = 0; t < T; t++)
				M_[i] = std::max( M_[i], y[ i * T + t] );
			assert( M_[i] < K_ );
		}
		//
		// set logfac_
		logfac_.resize(K_);
		logfac_[0]  = 0.0;
		logfac_[1]  = 0.0;
		for(size_t k = 2; k < K_; k++)
			logfac_[k] = log( double(k) ) + logfac_[k-1];
		//
		// set log_pik_
		log_pik_.resize(R_ * K_);
		for(size_t ell = 0; ell < R_ * K_; ell++)
			log_pik_[ell] = 0.0;
		implement_ran_likelihood(fixed_in, random_in, log_pik_);
	}
	// implementaion of ran_likelihood
	template <class Float>
	vector<Float> implement_ran_likelihood(
		const vector<Float>&  theta   ,
		const vector<Float>&  u       ,
		vector<Float>&        log_pik )
	{	assert( log_pik.size() == R_* K_ );
		vector<Float> vec(1);
		//
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
		// y_{i,t} * log( qt )
		vector<Float> yit_log_q(R_ * T_);
		// log( 1.0 - qt )
		vector<Float> log_1q(T_);
		for(size_t t = 0; t < T_; t++)
		{	Float ex  = exp( - u[t] - theta[1] );
			Float q   = one / (one + ex );
			log_1q[t] = log(one - q + eps);
			for(size_t i = 0; i < R_; i++)
				yit_log_q[ i * T_ + t ]  = Float(y_[i*T_+t]) * log(q + eps);
		}
		//
		// log( theta[0]^k * exp( - theta[0] ) / k! )
		vector<Float> log_poisson(K_);
		for(size_t k = 0; k < K_; k++)
		{	log_poisson[k] = log( theta[0] + eps ) * Float(k);
			log_poisson[k] -= theta[0];
			log_poisson[k] -= logfac_[k];
		}
		//
		// log[ p(y|theta, u) ] to vec[0]
		for(size_t i = 0; i < R_; i++)
		{	// compute maximum of log_pik for this i
			Float max_log_pik = log_pik[ i * K_ + M_[i] ];
			for(size_t k = M_[i] + 1; k < K_; k++)
			{	if( log_pik[ i * K_ + k ] > max_log_pik )
					max_log_pik = log_pik[ i * K_ + k ];
			}
			//
			// initialize p(y_i|theta,u)
			Float p_i = Float(0.0);
			for(size_t k = M_[i]; k < K_; k++)
			{	// compute p(y_i|N_i=k,theta,u) p(N_i=k|theta)
				//
				// initialize sum that needs to be calculated with Float
				// p(N_i=k|theta)
				Float float_sum = log_poisson[k];
				//
				// initilaize sum that does not need to use Float
				double double_sum = 0.0;
				//
				// now compute terms that depend on t
				for(size_t t = 0; t < T_; t++)
				{	size_t yit = y_[ i * T_ + t ];
					//
					// k choose yit = k! / [ yit! * (k - yit)! ]
					// log [ (k choose yit) ]
					double_sum += logfac_[k] - logfac_[yit] - logfac_[k - yit];
					//
					// log ( qt^yit )
					float_sum += yit_log_q[i * T_ + t];
					//
					// log [ (1 - qt)^(k - yit) ]
					float_sum += Float(k - yit) * log_1q[t];
				}
				// log[ p(y_i|N_i=k,theta,u) p(N_i=k|theta)
				// (output elements of log_pik are not affected by its input)
				log_pik[ i * K_ + k ] = float_sum + double_sum;
				//
				// normalization avoids exponential overflow
				// p_i += p(y_i|N_i=k,theta,u) p(N_i=k|theta) / exp(max_log_pik)
				p_i += exp( log_pik[ i * K_ + k ] - max_log_pik );
			}
			// vec[0] += log[ p(y_i|theta,u) ]
			vec[0] += log( p_i );
			// following term does not depend on the fixed or random effects
			// vec[0] += Float(K_ - M_[i]) * max_log_pik;
		}
		// - log [ p(y|theta,u) p(u|theta) ]
		vec[0] = - vec[0];
		//
		// result may be inifite or nan when only computing log_pik
		// (initial log_pik is used to compute normalization factors)
		return vec;
	}
// ------------------------------------------------------------------------
public:
	virtual vector<a2_double> ran_likelihood(
		const vector<a2_double>& fixed_vec  ,
		const vector<a2_double>& random_vec )
	{	vector<a2_double> log_pik(R_ * K_), vec(1);
		for(size_t ell = 0; ell < R_ * K_; ell++)
			log_pik[ell] = a2_double(log_pik_[ell]);
		vec = implement_ran_likelihood(fixed_vec, random_vec, log_pik);
		//
		// make sure result is finite
		double inf = std::numeric_limits<double>::infinity();
		assert( vec[0] < + a2_double(inf) );
		assert( vec[0] > - a2_double(inf) );
		//
		return vec;
	}
	//
	virtual vector<a1_double> ran_likelihood(
		const vector<a1_double>& fixed_vec  ,
		const vector<a1_double>& random_vec )
	{	vector<a1_double> log_pik(R_ * K_), vec(1);
		for(size_t ell = 0; ell < R_ * K_; ell++)
			log_pik[ell] = a1_double(log_pik_[ell]);
		vec = implement_ran_likelihood(fixed_vec, random_vec, log_pik);
		//
		// make sure result is finite
		double inf = std::numeric_limits<double>::infinity();
		assert( vec[0] < + a1_double(inf) );
		assert( vec[0] > - a1_double(inf) );
		//
		return vec;
	}
};
template <class Value>
void label_print(const char* label, const Value& value)
{	std::cout << std::setw(35) << std::left << label;
	std::cout << " = " << value << std::endl;
}
} // END_EMPTY_NAMESPACE

int main(int argc, char *argv[])
{	bool ok = true;
	using std::endl;
	using std::string;
	//
	const char* arg_name[] = {
		"random_seed",
		"number_fixed_samples",
		"number_locations",
		"number_times",
		"max_population",
		"mean_population",
		"mean_logit_probability",
		"std_logit_probability",
		"quasi_fixed",
		"random_constraint",
		"trace_optimize_fixed"
	};
	size_t n_arg = sizeof(arg_name)/sizeof(arg_name[0]);
	//
	if( size_t(argc) != 1 + n_arg )
	{	// print usage error message
		std::cerr << "usage: " << argv[0];
		for(size_t i = 0; i < n_arg; i++)
			std::cerr << " \\ \n\t" << arg_name[i];
		std::cerr << "\n";
		std::exit(1);
	}
	//
	// get command line arguments
	assert( n_arg == 11 );
	size_t random_seed            = std::atoi( argv[1] );
	size_t number_fixed_samples   = std::atoi( argv[2] );
	size_t number_locations       = std::atoi( argv[3] );
	size_t number_times           = std::atoi( argv[4] );
	size_t max_population         = std::atoi( argv[5] );
	double mean_population        = std::atof( argv[6] );
	double mean_logit_probability = std::atof( argv[7] );
	double std_logit_probability  = std::atof( argv[8] );
	//
	string quasi_fixed_str          = argv[9];
	string random_constraint_str    = argv[10];
	string trace_optimize_fixed_str = argv[11];
	bool quasi_fixed          =  quasi_fixed_str == "true";
	bool random_constraint    =  random_constraint_str == "true";
	bool trace_optimize_fixed =  trace_optimize_fixed_str == "true";
	//
	// print the commmand line with actual converted values
	std::cout   << argv[0]
	<< " " << random_seed
	<< " " << number_fixed_samples
	<< " " << number_locations
	<< " " << number_times
	<< " " << max_population
	<< " " << mean_population
	<< " " << mean_logit_probability
	<< " " << std_logit_probability
	<< " " << quasi_fixed_str
	<< " " << random_constraint_str
	<< " " << trace_optimize_fixed_str
	<< endl;
	//
	// print the command line arugments with labels for each value
	for(size_t i = 0; i < n_arg; i++)
		label_print(arg_name[i], argv[1+i]);
	//
	assert( number_fixed_samples > 0 );
	assert( number_locations > 0 );
	assert( number_times > 0 );
	assert( max_population > 0.0 );
	assert( mean_population > 0.0 );
	assert( std_logit_probability > 0.0 );
	assert( quasi_fixed || quasi_fixed_str=="false" );
	assert( random_constraint || random_constraint_str=="false" );
	assert( trace_optimize_fixed || trace_optimize_fixed_str=="false" );
	//
	// Get actual seed (different when random_seed is zero).
	size_t actual_seed     = CppAD::mixed::new_gsl_rng( random_seed );
	label_print("actual_seed", actual_seed);
	//
	// number of fixed and random effects
	size_t n_fixed  = 3;
	size_t n_random = number_times;
	//
	// number of random effects
	size_t T = number_times;
	size_t R = number_locations;
	size_t K = max_population;
	//
	vector<double> theta_sim(n_fixed);
	theta_sim[0] = mean_population;
	theta_sim[1] = mean_logit_probability;
	theta_sim[2] = std_logit_probability;
	//
	// simulate y
	vector<size_t> y(R * T);
	simulate(random_constraint, R, T, theta_sim, y);

	// lower and upper limits
	vector<double> fix_constraint_lower, fix_constraint_upper;
	vector<double>
		theta_lower(n_fixed), theta_in(n_fixed), theta_upper(n_fixed);
	//
	// theta[0] is estimate of mean_population
	theta_lower[0] = 0.0;
	theta_in[0]    = mean_population / 2.0;
	theta_upper[0] = 2.0 * mean_population;
	//
	// theta[1] is estimate of mean_logit_probability
	theta_lower[1] = mean_logit_probability - 2.0;
	theta_in[1]    = mean_logit_probability - 0.5;
	theta_upper[1] = mean_logit_probability + 2.0;
	//
	// theta[2] is estimate of std_logit_probability
	theta_lower[2] = std_logit_probability / 10.;
	theta_in[2]    = std_logit_probability / 2.0;
	theta_upper[2] = std_logit_probability * 10.;

	// random constraints
	CppAD::mixed::sparse_mat_info A_info;
	if( random_constraint )
	{	A_info.resize(T);
		for(size_t t = 0; t < T; t++)
		{	A_info.row[t] = 0;
			A_info.col[t] = t;
			A_info.val[t] = 1.0;
		}
	}

	// initialize random effects to start optimization at
	vector<double>  u_in(T);
	for(size_t t = 0; t < T; t++)
		u_in[t] = 0.0;
	//
	// start timing of initialization
	std::time_t start_time = std::time( CPPAD_MIXED_NULL_PTR );
	//
	// create derived object
	mixed_derived mixed_object(
		R, T, K, quasi_fixed, A_info, y, theta_in, u_in
	);
	//
	// initialize derived object
	std::map<string, size_t> size_map =
		mixed_object.initialize(theta_in, u_in);
	//
	// end timing of initilization
	std::time_t end_time = std::time( CPPAD_MIXED_NULL_PTR );
	label_print("initialization_seconds", end_time - start_time);
	//
	// print amoumt of memory added to mixed_object during initialize
	size_t num_bytes_before = size_map["num_bytes_before"];
	size_t num_bytes_after  = size_map["num_bytes_after"];
	string bytes_str = CppAD::to_string(num_bytes_after - num_bytes_before);
	size_t n_char = bytes_str.size();
	string initialization_bytes = "";
	for(size_t i = 0; i < n_char; i++)
	{	if( i != 0 && (n_char - i) % 3 == 0 )
			initialization_bytes += ",";
		initialization_bytes += bytes_str[i];
	}
	label_print("initialization_bytes", initialization_bytes);
	//
	// ipopt options for optimizing the random effects
	string random_ipopt_options =
		"String  sb                        yes\n"
		"String  derivative_test           none\n"
		"String  derivative_test_print_all no\n"
		"Numeric tol                       1e-8\n"
		"Integer max_iter                  40\n"
		"Integer print_level               0\n"
	;
	// ipopt options for optimizing the fixd effects
	string fixed_ipopt_options =
		"String  sb                        yes\n"
		"String  derivative_test           none\n"
		"String  derivative_test_print_all no\n"
		"Numeric tol                       1e-8\n"
		"Integer max_iter                  40\n"
	;
	if( trace_optimize_fixed )
		fixed_ipopt_options += "Integer print_level               5\n";
	else
		fixed_ipopt_options += "Integer print_level               0\n";
	;
	double inf = std::numeric_limits<double>::infinity();
	vector<double> u_lower(n_random), u_upper(n_random);
	for(size_t i = 0; i < n_random; i++)
	{	u_lower[i] = -inf;
		u_upper[i] = +inf;
	}
	//
	// optimize fixed effects
	start_time = std::time( CPPAD_MIXED_NULL_PTR );
	if( trace_optimize_fixed )
		std::cout << endl;
	CppAD::mixed::fixed_solution solution = mixed_object.optimize_fixed(
		fixed_ipopt_options,
		random_ipopt_options,
		theta_lower,
		theta_upper,
		fix_constraint_lower,
		fix_constraint_upper,
		theta_in,
		u_lower,
		u_upper,
		u_in
	);
	end_time = std::time( CPPAD_MIXED_NULL_PTR );
	if( trace_optimize_fixed )
		std::cout << endl;
	label_print("optimize_fixed_seconds", end_time - start_time);
	//
	start_time = std::time( CPPAD_MIXED_NULL_PTR );
	vector<double> theta_out = solution.fixed_opt;
	vector<double> u_out     = mixed_object.optimize_random(
		random_ipopt_options,
		theta_out,
		u_lower,
		u_upper,
		u_in
	);
	end_time = std::time( CPPAD_MIXED_NULL_PTR );
	label_print("optimize_random_seconds", end_time - start_time);
	//
	// sum of random effects
	double sum_random_effects = 0.0;
	for(size_t j = 0; j < n_random; j++)
		sum_random_effects += u_out[j];
	if( random_constraint )
		ok &= std::fabs( sum_random_effects ) < 1e-8;
	label_print("sum_random_effects", sum_random_effects);
	//
	// information matrix
	start_time = std::time( CPPAD_MIXED_NULL_PTR );
	CppAD::mixed::sparse_mat_info
	information_info = mixed_object.information_mat(solution, u_out);
	end_time = std::time( CPPAD_MIXED_NULL_PTR );
	label_print("information_mat_seconds", end_time - start_time);
	//
	start_time = std::time( CPPAD_MIXED_NULL_PTR );
	vector<double> sample( number_fixed_samples * n_fixed );
	mixed_object.sample_fixed(
		sample,
		information_info,
		solution,
		theta_lower,
		theta_upper,
		u_out
	);
	end_time = std::time( CPPAD_MIXED_NULL_PTR );
	label_print("sample_fixed_seconds", end_time - start_time);
	//
	// compute the sample standard deviations
	vector<double> sample_std(n_fixed), estimate_ratio(n_fixed);
	for(size_t j = 0; j < n_fixed; j++)
		sample_std[j] = 0.0;
	for(size_t i = 0; i < number_fixed_samples; i++)
	{	for(size_t j = 0; j < n_fixed; j++)
		{	double diff    = sample[i * n_fixed + j] - theta_out[j];
			sample_std[j] += diff * diff;
		}
	}
	for(size_t j = 0; j < n_fixed; j++)
	{	sample_std[j] = std::sqrt( sample_std[j] / number_fixed_samples );
		// check if results results are reasonable
		estimate_ratio[j] = ( theta_out[j] - theta_sim[j] ) / sample_std[j];
		ok  &= std::fabs(estimate_ratio[j]) < 4.0;
	}
	label_print("mean_population_estimate", theta_out[0]);
	label_print("mean_logit_probability_estimate", theta_out[1]);
	label_print("std_logit_probability_estimate", theta_out[2]);
	//
	label_print("mean_population_std", sample_std[0]);
	label_print("mean_logit_probability_std", sample_std[1]);
	label_print("std_logit_probability_std", sample_std[2]);
	//
	label_print("mean_population_ratio", estimate_ratio[0]);
	label_print("mean_logit_probability_ratio", estimate_ratio[1]);
	label_print("std_logit_probability_ratio", estimate_ratio[2]);
	//
	if( ok )
		label_print("capture_xam_ok", "true");
	else
		label_print("capture_xam_ok", "false");
	//
	// free memory allocated by new_gsl_rng
	CppAD::mixed::free_gsl_rng();
	if( ok )
		return 0;
	return 1;
}
// END C++
