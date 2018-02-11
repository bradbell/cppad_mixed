/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/ipopt_fixed.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
$begin ipopt_fixed_eval_grad_f$$
$spell
	CppAD
	ran_obj
	cppad
	obj
	ipopt
	bool
	eval
	const
$$

$section Compute Gradient of the Objective$$

$head Syntax$$
$icode%ok% = eval_grad_f(%n%, %x%, %new_x%, %grad_f%)%$$

$head n$$
is the number of variables in the problem (dimension of x).

$head x$$
is the value for the primal variables at which the gradient
$latex \nabla f(x)$$ is computed (has size $icode n$$).

$head new_x$$
if true, no Ipopt evaluation method was previous called with the same
value for $icode x$$.

$head grad_f$$
is set to the value for the gradient $latex \nabla f(x)$$
(has size $icode m$$).

$head ok$$
If set to false, the optimization will terminate with status set to
$cref/USER_REQUESTED_STOP
	/ipopt_fixed_finalize_solution/status/USER_REQUESTED_STOP/$$.

$head Prototype$$
$srccode%cpp% */
bool ipopt_fixed::eval_grad_f(
	Index           n         ,  // in
	const Number*   x         ,  // in
	bool            new_x     ,  // in
	Number*         grad_f    )  // out
/* %$$
$end
*/
{	try
	{	try_eval_grad_f(n, x, new_x, grad_f);
	for(size_t j = 0; j < size_t(n); j++)
		grad_f[j] *= scale_f_;
	}
	catch(const CppAD::mixed::exception& e)
	{	error_message_ = e.message("ipopt_fixed::eval_g");
		for(size_t j = 0; j < n_fixed_; j++)
			error_fixed_[j] = x[j];
		return false;
	}
	return true;
}
void ipopt_fixed::try_eval_grad_f(
	Index           n         ,  // in
	const Number*   x         ,  // in
	bool            new_x     ,  // in
	Number*         grad_f    )  // out
{
	assert( n > 0 && size_t(n) == n_fixed_ + fix_likelihood_nabs_ );
	//
	// fixed effects
	for(size_t j = 0; j < n_fixed_; j++)
		fixed_tmp_[j] = double( x[j] );
	//
	// random part of objective
	assert( H_beta_tmp_.size() == n_fixed_ );
	for(size_t j = 0; j < n_fixed_; j++)
		H_beta_tmp_[j] = Number(0.0);

	d_vector r_fixed(n_fixed_);
	if( n_random_ > 0 )
	{
		// compute the optimal random effects corresponding to fixed effects
		if( new_x )
			new_random(fixed_tmp_);
		// Jacobian for random part of the Lalpace objective
		mixed_object_.ran_obj_jac(
			fixed_tmp_, random_cur_, H_beta_tmp_
		);
	}
	//
	// Jacobian of fixed part of likelihood
	// (2DO: do not revaluate when eval_jac_g has same x)
	s_vector fix_like_jac_row = mixed_object_.fix_like_jac_.subset.row();
	s_vector fix_like_jac_col = mixed_object_.fix_like_jac_.subset.col();
	d_vector fix_like_jac_val = mixed_object_.fix_like_jac_.subset.val();
	mixed_object_.fix_like_jac(
		fixed_tmp_,
		fix_like_jac_row,
		fix_like_jac_col,
		fix_like_jac_val
	);
	//
	// Laplace objective part of grad_f
	for(size_t j = 0; j < n_fixed_; j++)
	{	assert( j < size_t(n) );
		grad_f[j] = Number( H_beta_tmp_[j] );
	}
	// auxillary variable part of grad_f
	for(size_t j = 0; j < fix_likelihood_nabs_; j++)
	{	assert( n_fixed_ + j < size_t(n) );
		grad_f[n_fixed_ + j] = Number( 1.0 );
	}
	// fixed likelihood part of grad_f
	for(size_t k = 0; k < fix_like_jac_row.size(); k++)
	{	if( fix_like_jac_row[k] == 0 )
		{	size_t j = fix_like_jac_col[k];
			assert( j < size_t(n) );
			grad_f[j] += Number( fix_like_jac_val[k] );
		}
	}
	for(size_t j = 0; j < size_t(n); j++)
	{	if( CppAD::isnan( grad_f[j] ) ) throw CppAD::mixed::exception(
			"", "gradient has a nan"
		);
	}
	//
	return;
}
} } // END_CPPAD_MIXED_NAMESPACE
