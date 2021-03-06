// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-20 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin sample_random$$
$spell
	vec
	gsl_rng
	cppad
	const
	CppAD
	std
	ipopt
$$

$section Simulation the Posterior Distribution for Random Effects$$

$head Syntax$$
$icode%error_msg% = %mixed_object%.sample_random(
	%sample%,
	%fixed_vec%,
	%random_ipopt_options%,
	%random_lower%,
	%random_upper%,
	%random_in%
)%$$

$head See Also$$
$cref sample_fixed$$

$head Prototype$$
$srcthisfile%0%// BEGIN PROTOTYPE%// END PROTOTYPE%1%$$

$head Public$$
This $code cppad_mixed$$ $cref base_class$$ member function is public.

$head Purpose$$
This routine draws samples from
the asymptotic posterior distribution for the
random effects given the model, the data, and the fixed effects; see
$cref/sparse observed information/theory/Sparse Observed Information/$$.

$head manage_gsl_rng$$
It is assumed that
$cref/get_gsl_rng/manage_gsl_rng/get_gsl_rng/$$ will return
a pointer to a GSL random number generator.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head sample$$
This argument has prototype
$codei%
	CppAD::vector<double>& %sample%
%$$
and its size is a multiple of
$cref/n_random/derived_ctor/n_random/$$.
The input value of its elements does not matter.
We define
$codei%
	%n_sample% = %sample_size% / %n_random%
%$$
If $icode error_msg$$ is empty, upon return
for $codei%i% = 0 , %...%, %n_sample%-1%$$,
$codei%j% = 0 , %...%, %n_random%-1%$$,
$codei%
	%sample%[ %i% * %n_random% + %j% ]
%$$
is the $th j$$ component of the $th i$$ sample of the
optimal random effects.
The statistics of these samples is specified under
$cref/covariance/sample_random/Covariance/$$ below.

$head random_ipopt_options$$
This argument has prototype
$codei%
	const std::string& %random_ipopt_options%
%$$
and is the $cref ipopt_options$$ for optimizing the random effects.

$head fixed_vec$$
This argument specifies the value of the
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$
vector $latex \theta$$.

$head random_lower$$
This argument must have size equal to
$cref/n_random/derived_ctor/n_random/$$ and
specifies the lower limits for the optimization of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$.
The value minus infinity can be used to specify no lower limit.

$head random_upper$$
This argument must have size equal to
$cref/n_random/derived_ctor/n_random/$$ and
specifies the upper limits for the optimization of the random effect.
The value plus infinity can be used to specify no lower limit.

$head random_in$$
This argument must have size equal to
$cref/n_random/derived_ctor/n_random/$$ and
specifies the initial value used for the optimization of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$.
It must hold that
$codei%
	%random_lower%[%i%] <= %random_in%[%i%] <= %random_upper%[%i%]
%$$
for each valid index $icode i$$.

$head Covariance$$
Each sample of the random effects is an independent normal.
The mean for this distribution is the
$cref/optimal random effects/theory/Optimal Random Effects, u^(theta)/$$
$latex \hat{u} ( \theta )$$.
The variance of this distribution
is the inverse of the observed information
matrix; i.e.
$latex \[
	f_{uu} [ \theta , \hat{u} ( \theta ) ] ^{-1}
\] $$
This normal distribution is censored to be within the limits
$icode random_lower$$, $icode random_upper$$.

$head error_msg$$
If $icode error_msg$$ is empty (non-empty),
$cref/sample/sample_random/sample/$$
values have been calculated (have not been calculated).
If $icode error_msg$$ is non-empty,
it is a message describing the problem.


$children%example/user/sample_random.cpp
%$$
$head Example$$
The file $cref sample_random.cpp$$ is an example
and test of $code sample_random$$.
$end
-----------------------------------------------------------------------------
*/

# include <Eigen/Core>
# include <Eigen/Cholesky>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/manage_gsl_rng.hpp>
# include <gsl/gsl_randist.h>
# include <cppad/mixed/exception.hpp>

std::string cppad_mixed::try_sample_random(
	d_vector&          sample               ,
	const std::string& random_ipopt_options ,
	const d_vector&    fixed_vec            ,
	const d_vector&    random_lower         ,
	const d_vector&    random_upper         ,
	const d_vector&    random_in            )
{	// case where there is nothing to do
	if( n_random_ == 0 )
		return "";
	//
	assert( sample.size() % n_random_ == 0   );
	assert( fixed_vec.size()    == n_fixed_  );
	assert( random_lower.size() == n_random_ );
	assert( random_upper.size() == n_random_ );
	assert( random_in.size()    == n_random_ );
	//
	// number of samples
	size_t n_sample = sample.size() / n_random_;
	//
	// optimal random effects
	d_vector random_opt;
	random_opt = optimize_random(
		random_ipopt_options, fixed_vec, random_lower, random_upper, random_in
	);
	// update the Cholesky factor corresponding to f_uu (theta, u)
	update_factor(fixed_vec, random_opt);
	//
	for(size_t i_sample = 0; i_sample < n_sample; i_sample++)
	{	// simulate a normal with mean zero and variance one
		d_vector w(n_random_);
		for(size_t j = 0; j < n_random_; j++)
			w[j] = gsl_ran_gaussian(CppAD::mixed::get_gsl_rng(), 1.0);
		//
		// set v to cholesky factor of f_uu(theta, u)^{-1} times w
		d_vector v(n_random_);
		bool ok = ldlt_ran_hes_.sim_cov(w, v);
		if( ! ok )
		{	std::string msg = "sample_random: Hessian w.r.t random effects"
				" is not positive definite";
			return msg;
		}
		//
		// add random_opt an truncate to random limits
		for(size_t j = 0; j < n_random_; j++)
		{	double samp = random_opt[j] + v[j];
			samp = std::min(samp, random_upper[j]);
			samp = std::max(samp, random_lower[j]);
			sample[i_sample * n_random_ + j] = samp;
		}
	}
	return "";
}
// --------------------------------------------------------------------------
// BEGIN PROTOTYPE
std::string cppad_mixed::sample_random(
	d_vector&          sample               ,
	const std::string& random_ipopt_options ,
	const d_vector&    fixed_vec            ,
	const d_vector&    random_lower         ,
	const d_vector&    random_upper         ,
	const d_vector&    random_in            )
// END PROTOTYPE
{	std::string error_msg = "";
	try
	{	error_msg = try_sample_random(
			sample                ,
			random_ipopt_options  ,
			fixed_vec             ,
			random_lower          ,
			random_upper          ,
			random_in
		);
	}
	catch(const CppAD::mixed::exception& e)
	{	error_msg = e.message("sample_random");
	}
	return error_msg;
}
