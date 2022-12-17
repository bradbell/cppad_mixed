/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/ipopt_fixed.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
$begin ipopt_fixed_get_starting_point$$
$spell
	CppAD
	ran_obj
	cppad
	obj
	init
	ipopt
	bool
$$

$nospell
$bold This is old cppad_mixed documentation:$$ Here is a link to its
$href%http://bradbell.github.io/cppad_mixed%current documentation%$$.
$$
$section Return Initial Values Where Optimization is Started$$

$head Syntax$$
$icode%ok% = get_starting_point(
	%n%, %init_x%, %x%, %init_z%, %z_L%, %z_U%, %m%, %init_lambda%, %lambda%
)%$$

$head n$$
is the number of variables in the problem (dimension of x).

$head init_x$$
assumed true which means the ipopt options specify that the this routine
will provide an initial value for $icode x$$.

$head x$$
if $icode init_x$$ is true,
set to the initial value for the primal variables (has size $icode n$$).

$head init_z$$
assumes $icode init_z$$ is false.
If it were true, the ipopt options specify that the this routine
will provide an initial value for $icode x$$ upper and lower bound
multipliers.

$head z_L$$
if $icode init_z$$ is true,
set to the initial value for the lower bound multipliers (has size $icode n$$).

$head z_U$$
if $icode init_z$$ is true,
set to the initial value for the upper bound multipliers (has size $icode n$$).

$head init_lambda$$
assumes $icode init_lambda$$ is false.
If it were true, the ipopt options specify that the this routine
will provide an initial value for $icode g(x)$$ upper and lower bound
multipliers.

$head lambda$$
if $icode init_lambda$$ is true,
set to the initial value for the $icode g(x)$$ multipliers
(has size $icode m$$).

$head ok$$
If set to false, the optimization will terminate with status set to
$cref/USER_REQUESTED_STOP
	/ipopt_fixed_finalize_solution/status/USER_REQUESTED_STOP/$$.

$head Prototype$$
$srccode%cpp% */
bool ipopt_fixed::get_starting_point(
	Index           n            ,  // in
	bool            init_x       ,  // in
	Number*         x            ,  // out
	bool            init_z       ,  // in
	Number*         z_L          ,  // out
	Number*         z_U          ,  // out
	Index           m            ,  // in
	bool            init_lambda  ,  // in
	Number*         lambda       )  // out
/* %$$
$end
*/
{	assert( adaptive_called_ == true );
	assert( init_x == true );
	assert( n > 0 );
	assert( size_t(n) == n_fixed_ + fix_likelihood_nabs_ );
	assert( m >= 0 );
	assert( size_t(m) == 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_ );
	//
	if( warm_start_.x_info.size() == 0 )
	{	assert( init_z == false );
		assert( init_lambda == false );

		// use input values for fixed effects
		for(size_t j = 0; j < n_fixed_; j++)
			x[j] = fixed_in_[j];

		// set auxillary variables to corresponding minimum feasible value
		for(size_t j = 0; j < fix_likelihood_nabs_; j++)
			x[n_fixed_ + j] = std::fabs( fix_likelihood_vec_tmp_[1 + j] );
	}
	else
	{	assert( warm_start_.x_info.size() == size_t(n) );
		assert( warm_start_.g_info.size() == size_t(m) );
		assert( init_z  == true );
		assert( init_lambda == true );
		for(Index j = 0; j < n; ++j)
		{	x[j]   = warm_start_.x_info[j].x;
			z_L[j] = warm_start_.x_info[j].z_L;
			z_U[j] = warm_start_.x_info[j].z_U;
		}
		for(Index i = 0; i < m; ++i)
			lambda[i] = warm_start_.g_info[i].lambda;
	}
	return true;
}
} } // END_CPPAD_MIXED_NAMESPACE
