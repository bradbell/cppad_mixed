/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/configure.hpp>
# include <cppad/mixed/ipopt_fixed.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
$begin ipopt_fixed_get_scaling_parameters$$
$spell
	Ipopt
$$

$section Inform Ipopt of Scaling Parameters$$

This is not documented.

$end
*/
bool ipopt_fixed::get_scaling_parameters(
	Number&            obj_scaling    ,
	bool&              use_x_scaling  ,
	Index              n              ,
	Number*            x_scaling      ,
	bool&              use_g_scaling  ,
	Index              m              ,
	Number*            g_scaling      )
{
# if CPPAD_MIXED_HIDE_IPOPT_SCALING
	assert( false );
	return false;
# else
	assert( n > 0 );
	assert( size_t(n) == n_fixed_ + fix_likelihood_nabs_ );
	assert( m >= 0 );
	assert( size_t(m) == 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_ );
	//
	obj_scaling   = Number(scale_f_);
	//
	use_x_scaling = false;
	for(int j = 0; j < n; j++)
		x_scaling[j] = Number(1.0);
	//
	use_g_scaling = true;
	for(int i = 0; i < m; i++)
		g_scaling[i] = Number( scale_g_[i] );
	//
	return true;
# endif
}
} } // END_CPPAD_MIXED_NAMESPACE
