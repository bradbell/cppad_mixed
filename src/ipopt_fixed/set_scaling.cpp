/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/ipopt_fixed.hpp>
# include <cppad/mixed/one_dim_derivative_chk.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

/*
$begin ipopt_fixed_set_scaling$$
$spell
$$

$section Set Scaling Factors$$

$head Syntax$$
$icode%ok% = set_scaling(
	%x_scale%, %x_lower%, %x_upper%, %grad_f%, %jac_g%)

$head x_scale$$
is the argument value at which the scaling is computer.

$head x_lower$$
is the optimization lower limit for x.

$head x_upper$$
is the optimization upper limit for x.

$head grad_f$$
is the gradient of f(x) at $icode x_scale$$.

$head jac_g$$
is the Jacobian of g(x) at $icode x_scale$$.
This is a sparse representation of the jacobian using the member variables
$code jac_g_row_$$ and $code jac_g_col_$$.

$head scale_f_$$
The input value of this member varable must be one.
Upon return it is set to a value that should be used to multiply
the values return by $code eval_f$$.

$head scale_g_$$
The input value of this member varable must be a vector with size
equal to the range space of g(x) and all its components must be one.
Upon return $codei%scale_g_[%i%]%$$ is set to a value that should be used to
multiply the $th i$$ component returned by $code eval_g$$.

$head Prototype$$
$srcthisfile%0%// BEGIN_PROTOTYPE%// END_PROTOTYPE%1%$$

$end
-------------------------------------------------------------------------------
*/
// BEGIN_PROTOTYPE
bool ipopt_fixed::set_scaling(
	const d_vector& x_scale ,
	const d_vector& x_lower ,
	const d_vector& x_upper ,
	const d_vector& grad_f  ,
	const d_vector& jac_g   )
// END_PROTOTYPE
{	//
	// m: number of components if g
	size_t m    = 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_;
	//
# ifndef NDEBUG
	assert( scale_f_ == 1.0 );
	assert( scale_g_.size() == m );
	for(size_t i = 0; i < m; ++i)
		assert( scale_g_[i] == 1.0 );
# endif
	//
	// scale_max: maximum scaling factor
	const double scale_max = 1e+14;
	//
	// scale_min: minimum scaling factor
	double scale_min       = 1.0 / scale_max;
	//
	// max_jac_g : maximum absolute partial for each component of g
	// not counting components with equal lower and upper limits
	d_vector max_jac_g(m);
	for(size_t i = 0; i < m; i++)
		max_jac_g[i] = scale_min;
	for(size_t k = 0; k < nnz_jac_g_; k++)
	{	size_t i = jac_g_row_[k];
		size_t j = jac_g_col_[k];
		if( x_lower[j] < x_upper[j] )
			max_jac_g[i] = std::max( max_jac_g[i], std::fabs( jac_g[k] ) );
	}
	//
	// max_grad_f : maximum absolute partial for objective.
	// Skip absolute value partials in f.
	double max_grad_f = scale_min;
	// fixed effect terms
	for(size_t j = 0; j < n_fixed_; j++) if( x_lower[j] < x_upper[j] )
		max_grad_f = std::max( max_grad_f, std::fabs( grad_f[j] ) );
	// skip axuillary variable terms (gradients are one)
# ifndef NDEBUG
	for(size_t j = 0; j < fix_likelihood_nabs_; j++)
		assert( grad_f[n_fixed_ + j] == 1.0 );
# endif
	// include absolute value terms in g(x) in max_grad_f
	for(size_t i = 0; i < fix_likelihood_nabs_; i++)
	{	assert( max_jac_g[2*i] == max_jac_g[2*i+1] );
		max_grad_f = std::max( max_grad_f, std::fabs( max_jac_g[2*i] ) );
	}
	//
	// set scale_f_
	scale_f_ = std::max( scale_min, 1.0 / max_grad_f);
	scale_f_ = std::min( scale_max, scale_f_);
	//
	// set scale_g_
	assert( scale_g_.size() == m );
	for(size_t i = 0; i < m; i++)
	{	if( i < 2 * fix_likelihood_nabs_ )
			scale_g_[i] = scale_f_;
		else
		{	scale_g_[i] = std::max( scale_min, 1.0 / max_jac_g[i] );
			scale_g_[i] = std::min( scale_max, scale_g_[i] );
		}
	}
	//
	return true;
}
} } // END_CPPAD_MIXED_NAMESPACE
