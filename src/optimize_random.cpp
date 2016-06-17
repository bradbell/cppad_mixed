// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/ipopt/solve.hpp>
# include <cppad/mixed/cppad_mixed.hpp>
/*
$begin optimize_random$$
$spell
	CppAD
	cppad
	vec
	const
	CppAD
	std
	Ipopt
	optimizer
$$

$section Optimize Random Effects$$

$head Syntax$$
$icode%random_out% =%$$
$icode%mixed_object%.optimize_random(
	%options%, %fixed_vec%, %random_lower%, %random_upper%, %random_in%
)%$$

$head Public$$
This $code cppad_mixed$$ member function is $cref public$$.

$head Purpose$$
This routine maximizes the
$cref/random likelihood/ran_likelihood/$$
corresponding to the object $icode mixed_object$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head options$$
This argument  has the prototype
$codei%
	const std::string& %options%
%$$
This is the $cref ipopt_options$$ for optimizing the random effects.

$head fixed_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fixed_vec%
%$$
It specifies the value of the
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$
vector $latex \theta$$.

$head random_lower$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_lower%
%$$
It must have size equal to
$cref/n_random/derived_ctor/n_random/$$ and
specifies the lower limits for the optimization of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$.
The value minus infinity can be used to specify no lower limit.

$head random_upper$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_upper%
%$$
It must have size equal to
$cref/n_random/derived_ctor/n_random/$$ and
specifies the upper limits for the optimization of the random effect.
The value plus infinity can be used to specify no lower limit.

$head random_in$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_in%
%$$
It must have size equal to
$cref/n_random/derived_ctor/n_random/$$ and
specifies the initial value used for the optimization of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$ vector $latex u$$.
It must hold that
$codei%
	%random_lower%[%i%] <= %random_in%[%i%] <= %random_upper%[%i%]
%$$
for each valid index $icode i$$.

$head random_out$$
The return value  has prototype
$codei%
	CppAD::vector<double> %random_out%
%$$
It is the final value (obtained by optimization) of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$.

$children%
	example/user/optimize_random_xam.cpp
%$$
$head Example$$
The file $cref optimize_random_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
// ----------------------------------------------------------------------------
// helper class used by optimize_random
// (not in empty namespace so can be friend of cppad_mixed class)
class optimize_random_ipopt {
public:
	// same as cppad_mixed::d_vector
	typedef CppAD::vector<double>              Dvector;

	// same as cppad_mixed::a1d_vector
	typedef CppAD::vector< CppAD::AD<double> > ADvector;
private:
	const size_t     n_abs_;
	const size_t     n_fixed_;
	const size_t     n_random_;
	ADvector         fixed_vec_;
	cppad_mixed&    mixed_object_;
public:
	// constructor
	optimize_random_ipopt(
		size_t           n_abs           ,
		size_t           n_random        ,
		const Dvector&   fixed_vec       ,
		cppad_mixed&    mixed_object
	) :
	n_abs_        ( n_abs )            ,
	n_fixed_      ( fixed_vec.size() ) ,
	n_random_     ( n_random )         ,
	mixed_object_( mixed_object )
	{	fixed_vec_.resize( n_fixed_ );
		for(size_t i = 0; i < n_fixed_; i++)
			fixed_vec_[i] = fixed_vec[i];
	}

	// evaluate objective and constraints
	void operator()(ADvector& fg, const ADvector& x)
	{	assert( fg.size() == 1 + 2 * n_abs_ );

		// extract the random effects from x
		ADvector random_vec(n_random_);
		for(size_t j = 0; j < n_random_; j++)
			random_vec[j] = x[j];

		// compute log-density vector
		ADvector both_vec(n_fixed_ + n_random_);
		mixed_object_.pack(fixed_vec_, random_vec, both_vec);
		ADvector vec = mixed_object_.ran_like_a1fun_.Forward(0, both_vec);

		// initialize smooth part of negative log-likelihood
		size_t k = 0;
		fg[k++]  = vec[0];

		// terms corresponding to data likelihood absolute values
		size_t n_abs   = vec.size() - 1;
		for(size_t j = 0; j < n_abs; j++)
		{	// x[ n_random_ + j] >= abs(log_den[1 + j]
			fg[k++] = x[ n_random_ + j] - vec[1 + j];
			fg[k++] = x[ n_random_ + j] + vec[1 + j];
			//
			// smooth contribution to log-likelihood
			fg[0]  += x[ n_random_ + j];
		}
	}
};
} } // END_CPPAD_MIXED_NAMESPACE

// ----------------------------------------------------------------------------
// optimize_random
CppAD::vector<double> cppad_mixed::optimize_random(
	const std::string& options         ,
	const d_vector&    fixed_vec       ,
	const d_vector&    random_lower    ,
	const d_vector&    random_upper    ,
	const d_vector&    random_in       )
{
	// make sure initialize has been called
	if( ! initialize_done_ )
	{	std::string error_message =
		"cppad_mixed::initialize was not called before optimize_random";
		fatal_error(error_message);
	}
	if( ! init_ran_like_done_ )
	{	std::string error_message =
		"cppad_mixed::optimize_random there are no random effects";
		fatal_error(error_message);
	}

	// number of fixed and random effects
	assert( n_fixed_  == fixed_vec.size() );
	assert( n_random_ == random_in.size() );
	assert( n_random_ == random_lower.size() );
	assert( n_random_ == random_upper.size() );

	// infinity
	double inf = std::numeric_limits<double>::infinity();

	// determine initial density vector
	d_vector both_vec(n_fixed_ + n_random_);
	pack(fixed_vec, random_in, both_vec);
	d_vector vec = ran_like_fun_.Forward(0, both_vec);

	// number of absolute value terms in objective
	size_t n_abs = vec.size() - 1;

	// number of independent variable is number of random effects
	// plus number of log-density terms that require absolute values
	size_t nx = n_random_ + n_abs;

	// set initial x vector
	d_vector xi(nx);
	for(size_t j = 0; j < n_random_; j++)
		xi[j] = random_in[j];
	for(size_t j = 0; j < n_abs; j++)
		xi[n_random_ + j] = CppAD::abs( vec[j] );

	// ipopts default value for infinity (use options to change it)
	double ipopt_infinity = 1e19;

	// set lower and upper limits for x
	d_vector xl(nx), xu(nx);
	for(size_t j = 0; j < nx; j++)
	{	if( random_lower[j] == - inf )
			xl[j] = - ipopt_infinity;
		else
			xl[j] = random_lower[j];
		if( random_lower[j] == + inf )
			xu[j] = + ipopt_infinity;
		else
			xu[j] = random_upper[j];
	}

	// set limits for abs contraint representation
	d_vector gl(2*n_abs), gu(2*n_abs);
	for(size_t j = 0; j < 2*n_abs; j++)
	{	gl[j] = 0.0;
		gu[j] = + ipopt_infinity;
	}

	// construct fg_eval  object
	CppAD::mixed::optimize_random_ipopt fg_eval(
		n_abs, n_random_, fixed_vec, *this
	);

	// optimizer options
	assert( options[ options.size() - 1] == '\n' );
	std::string solve_options = options + "Sparse  true  reverse \n";

	// return solution
	CppAD::ipopt::solve_result<d_vector> solution;

	// solve the optimization problem
	CppAD::ipopt::solve(
		solve_options, xi, xl, xu, gl, gu, fg_eval, solution
	);

	return solution.x;
}
