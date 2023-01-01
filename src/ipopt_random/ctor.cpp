/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/ipopt_random.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
$begin ipopt_random_ctor$$
$spell
	Ipopt
	CppAD
	vec
	cppad
	nlp
	inf
	nnz
$$

$nospell
$bold This is cppad_mixed--20220519 documentation:$$ Here is a link to its
$href%https://cppad-mixed.readthedocs.io%current documentation%$$.
$$
$section Ipopt Random Optimization Callback Constructor$$

$head Syntax$$
$codei%CppAD::mixed::ipopt_random %ipopt_object%(
	%fixed_vec%,
	%random_lower%,
	%random_upper%,
	%random_in%,
	%mixed_object%
)%$$

$head fixed_vec$$
specifies the value of the
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$
vector $latex \theta$$.
It is stored as a reference so it must exist for as long as
$icode ipopt_object$$ exists.

$head random_lower$$
this vector has size
$cref/n_random/derived_ctor/n_random/$$
and specifies the lower limits for the optimization of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$.
It is stored as a reference so it must exist for as long as
$icode ipopt_object$$ exists.
The value minus infinity can be used to specify no lower limit.

$head random_upper$$
this vector has size $icode n_random$$
and specifies the upper limits for the optimization of the random effects.
It is stored as a reference so it must exist for as long as
$icode ipopt_object$$ exists.
The value plus infinity can be used to specify no upper limit.

$head random_in$$
this vector has size $icode n_random$$ and specifies
the initial value used for the optimization of the random effects.
It is stored as a reference so it must exist for as long as
$icode ipopt_object$$ exists.
The value plus infinity can be used to specify no upper limit.

$head mixed_object$$
The argument $icode mixed_object$$ is an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head Member Variables$$

$subhead nlp_lower_bound_inf_$$
set to a finite value that is used by Ipopt for minus infinity.

$subhead nlp_upper_bound_inf_$$
set to a finite value that is used by Ipopt for plus infinity.

$subhead nnz_h_lag_$$
set to the number of non-zero entries in the Hessian of the Lagrangian.
This is the same as for the Hessian of the objective because there
are no constraints (except for box constraints) in this problem.

$srccode%cpp% */
ipopt_random::ipopt_random(
	const d_vector&     fixed_vec          ,
	const d_vector&     random_lower       ,
	const d_vector&     random_upper       ,
	const d_vector&     random_in          ,
	cppad_mixed&        mixed_object       ) :
/* %$$
$end
*/
n_fixed_             ( fixed_vec.size()  )                ,
n_random_            ( random_lower.size() )              ,
fixed_vec_           ( fixed_vec )                        ,
random_lower_        ( random_lower )                     ,
random_upper_        ( random_upper )                     ,
random_in_           ( random_in )                        ,
nnz_h_lag_           ( mixed_object.ran_hes_uu_rcv_.nnz() )  ,
mixed_object_        ( mixed_object    )                  ,
error_random_        ( n_random_ )
{	// -----------------------------------------------------------------------
	// set nlp_lower_bound_inf_, nlp_upper_bound_inf_
	double inf           = std::numeric_limits<double>::infinity();
	nlp_lower_bound_inf_ = - 1e19;
	nlp_upper_bound_inf_ = + 1e19;
	for(size_t j = 0; j < n_random_; j++)
	{	if( random_lower[j] != - inf ) nlp_lower_bound_inf_ =
				std::min(nlp_lower_bound_inf_, 1.1 * random_lower[j] );
		//
		if( random_upper[j] != inf ) nlp_upper_bound_inf_ =
				std::max(nlp_upper_bound_inf_, 1.1 * random_upper[j] );
	}
	objective_current_ = std::numeric_limits<double>::quiet_NaN();
	return;
}
} } // END_CPPAD_MIXED_NAMESPACE
