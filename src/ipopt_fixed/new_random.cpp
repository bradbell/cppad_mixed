/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/ipopt_fixed.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
$begin ipopt_fixed_new_random$$
$spell
	vec
	ipopt
	const
$$

$nospell
$bold This is cppad_mixed--20220519 documentation:$$ Here is a link to its
$href%https://cppad-mixed.readthedocs.io%current documentation%$$.
$$
$section Compute New Random Effects and Update Factor$$

$head Syntax$$
$codei%new_random(%fixed_vec%)%$$

$head ipopt_fixed$$
This is a private member function of the $cref ipopt_fixed$$ class.

$head n_random_$$
Is assumed that this member variable is greater than zero.

$head random_ipopt_options_$$
This member variable contains
the value of the
$cref/random_ipopt_options/ipopt_fixed_ctor/random_ipopt_options/$$
in the $code ipopt_fixed$$ constructor.

$head random_lower_$$
This member variable contains
the value of the $cref/random_lower/ipopt_fixed_ctor/random_lower/$$
in the $code ipopt_fixed$$ constructor.

$head random_upper_$$
This member variable contains
the value of the $cref/random_upper/ipopt_fixed_ctor/random_upper/$$
in the $code ipopt_fixed$$ constructor.

$head random_in_$$
This member variable contains
the value of the $cref/random_in/ipopt_fixed_ctor/random_in/$$
in the $code ipopt_fixed$$ constructor.

$head fixed_vec$$
This argument has prototype
$codei%
	const d_vector& %fixed_vec%
%$$
it is the value of the fixed effects that we are computing the random
effects and updated factor for.

$head random_cur_$$
This member variable contain is set the optimal random effects
corresponding to $icode fixed_vec$$.

$head mixed_object_$$
The factor in this member variables is updated using the call
$codei%
	mixed_object_.update_factor(%fixed_vec%, random_cur_)
%$$
see $cref update_factor$$ for side effects.

$head Prototype$$
$srccode%cpp% */
void ipopt_fixed::new_random(const d_vector& fixed_vec)
/* %$$
$end
*/
{	assert( n_random_ > 0 );
	// Compute the optimal random effects corresponding to fixed effects.
	// Use try_optimize_random instead of optimize_random, so that thows
	// are caught at the fixed effects, not random effects, level.
	random_cur_ = mixed_object_.try_optimize_random(
		random_ipopt_options_,
		fixed_vec,
		random_lower_,
		random_upper_,
		random_in_
	);
	mixed_object_.update_factor(fixed_vec, random_cur_);
}
} } // END_CPPAD_MIXED_NAMESPACE
