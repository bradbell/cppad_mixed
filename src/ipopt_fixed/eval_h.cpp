/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/configure.hpp>
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/ipopt_fixed.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
$begin ipopt_fixed_eval_h$$
$spell
	CppAD
	ran_obj
	cppad
	obj
	ipopt
	bool
	eval
	const
	obj
	nele_hess
	nnz
$$

$section Compute the Hessian of the Lagrangian$$

$head Syntax$$
$icode%ok% = eval_h(
	%n%, %x%, %new_x%,%obj_factor%, %m%, %lambda%, %new_lambda%,%$$
$icode%nele_hess%, %iRow%, %jCol%, %values%
)%$$

$head Lagrangian$$
The Lagrangian is defined to be
$latex \[
	L(x) = \alpha f(x) + \sum_{i=0}^{m-1} \lambda_i g_i (x)
\] $$

$head mixed_object.quasi_fixed_$$
It is assumed that this member variable is false.

$head n$$
is the number of variables in the problem (dimension of x).

$head x$$
is the value for the primal variables at which the
Hessian of the Lagrangian is computed (has size $icode n$$).

$head new_x$$
if true, no Ipopt evaluation method was previous called with the same
value for $icode x$$.

$head obj_factor$$
is the factor $latex \alpha$$ that multiplies the objective f(x)
in the definition of the Lagrangian.

$head m$$
is the number of constraints in the problem (dimension of g(x)).

$head lambda$$
is the value of the constraint multipliers $latex \lambda$$
at which the Hessian is to be evaluated (has size $icode m$$).

$head new_lambda$$
if true, no Ipopt evaluation method was previous called with the same
value for $icode lambda$$.

$head nele_hess$$
is the number of non-zero elements in the Hessian $latex L_{x,x} (x)$$; i.e.,
the same as
$cref/nnz_h_lag/ipopt_fixed_get_nlp_info/nnz_h_lag/$$.

$head iRow$$
If $icode values$$ is $code NULL$$,
$icode iRow$$ has size $icode nele_hess$$ and is set to the
row indices for the non-zero entries in the
lower triangle of the Hessian $latex L_{x,x} (x)$$.

$head jCol$$
If $icode values$$ is $code NULL$$,
$icode jCol$$ has size $icode nele_hess$$ and is set to the
column indices for the non-zero entries in the
lower triangle of the Hessian $latex L_{x,x} (x)$$.

$head values$$
If $icode values$$ is not $code NULL$$,
it has size $icode nele_hess$$ and $icode%values%[%k%]%$$
is set to the value of element of the Hessian $latex L_{x,x} (x)$$
with row index $icode%iRow%[%k%]%$$
and column index $icode%jCol%[%k%]%$$.

$head ok$$
if set to false, the optimization will terminate with status set to
$cref/USER_REQUESTED_STOP
	/ipopt_fixed_finalize_solution/status/USER_REQUESTED_STOP/$$.

$head Prototype$$
$srccode%cpp% */
bool ipopt_fixed::eval_h(
	Index         n              ,  // in
	const Number* x              ,  // in
	bool          new_x          ,  // in
	Number        obj_factor     ,  // in
	Index         m              ,  // in
	const Number* lambda         ,  // in
	bool          new_lambda     ,  // in
	Index         nele_hess      ,  // in
	Index*        iRow           ,  // out
	Index*        jCol           ,  // out
	Number*       values         )  // out
/* %$$
$end
*/
{	try
	{	try_eval_h(
			n,
			x,
			new_x,
			obj_factor,
			m,
			lambda,
			new_lambda,
			nele_hess,
			iRow,
			jCol,
			values
		);
	}
	catch(const CppAD::mixed::exception& e)
	{	error_message_ = e.message("ipopt_fixed::eval_h");
		for(size_t j = 0; j < n_fixed_; j++)
			error_fixed_[j] = x[j];
		return false;
	}
	return true;
}
void ipopt_fixed::try_eval_h(
	Index         n              ,  // in
	const Number* x              ,  // in
	bool          new_x          ,  // in
	Number        obj_factor     ,  // in
	Index         m              ,  // in
	const Number* lambda         ,  // in
	bool          new_lambda     ,  // in
	Index         nele_hess      ,  // in
	Index*        iRow           ,  // out
	Index*        jCol           ,  // out
	Number*       values         )  // out
{
	assert( ! mixed_object_.quasi_fixed_ );
	assert( n > 0 );
	assert( size_t(n) == n_fixed_ + fix_likelihood_nabs_ );
	assert( m >= 0 );
	assert( size_t(m) == 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_ );
	assert( size_t(nele_hess) == nnz_h_lag_ );
	if( values == NULL )
	{	for(size_t k = 0; k < nnz_h_lag_; k++)
		{	iRow[k] = Index( lag_hes_row_[k] );
			jCol[k] = Index( lag_hes_col_[k] );
		}
		return;
	}
	//
	// fixed effects
	for(size_t j = 0; j < n_fixed_; j++)
		fixed_tmp_[j] = double( x[j] );
	//
	// initialize return value
	for(size_t k = 0; k < nnz_h_lag_; k++)
		values[k] = Number( 0.0 );
	//
	// random part of objective
	if( n_random_ > 0 )
	{
		// compute the optimal random effects corresponding to fixed effects
		if( new_x )
			new_random(fixed_tmp_);
		// compute Hessian of random part of objective w.r.t. fixed effects
		w_laplace_obj_tmp_[0] = obj_factor;
# if CPPAD_MIXED_HIDE_IPOPT_SCALING
		w_laplace_obj_tmp_[0] = scale_f_ * obj_factor;
# endif
		//
		// include random constraints in this Hessian calculation
		size_t offset = 2 * fix_likelihood_nabs_ + n_fix_con_;
		for(size_t i = 0; i < n_ran_con_; i++)
		{	w_laplace_obj_tmp_[i+1] = lambda[offset + i];
# if CPPAD_MIXED_HIDE_IPOPT_SCALING
			w_laplace_obj_tmp_[i+1] = scale_g_[offset+i] * lambda[offset+i];
# endif
		}
		//
		mixed_object_.laplace_obj_hes(
			fixed_tmp_,
			random_cur_,
			w_laplace_obj_tmp_,
			laplace_obj_hes_info_.row,
			laplace_obj_hes_info_.col,
			laplace_obj_hes_info_.val
		);
		for(size_t k = 0; k < laplace_obj_hes_info_.row.size(); k++)
		{	size_t index = laplace_obj_hes_2_lag_[k];
			assert( index < nnz_h_lag_ );
			values[index] += Number( laplace_obj_hes_info_.val[k] );
		}
	}
	//
	// Hessian of Lagrangian of weighted fixed likelihood
	w_fix_likelihood_tmp_[0] = obj_factor;
# if CPPAD_MIXED_HIDE_IPOPT_SCALING
	w_fix_likelihood_tmp_[0] = scale_f_ * obj_factor;
# endif
	//
	for(size_t j = 0; j < fix_likelihood_nabs_; j++)
	{
		w_fix_likelihood_tmp_[1 + j] = lambda[2*j + 1] - lambda[2*j];
# if CPPAD_MIXED_HIDE_IPOPT_SCALING
		w_fix_likelihood_tmp_[1 + j] = scale_g_[2*j + 1] * lambda[2*j + 1]
		                             - scale_g_[2*j]     * lambda[2*j];
# endif
	}
	s_vector fix_like_hes_row = mixed_object_.fix_like_hes_.subset.row();
	s_vector fix_like_hes_col = mixed_object_.fix_like_hes_.subset.col();
	d_vector fix_like_hes_val = mixed_object_.fix_like_hes_.subset.val();
	mixed_object_.fix_like_hes(
		fixed_tmp_,
		w_fix_likelihood_tmp_,
		fix_like_hes_row,
		fix_like_hes_col,
		fix_like_hes_val
	);
	for(size_t k = 0; k < fix_like_hes_row.size(); k++)
	{	size_t index = fix_like_hes_2_lag_[k];
		assert( index < nnz_h_lag_ );
		values[index] += Number( fix_like_hes_val[k] );
	}
	//
	// Hessian of Lagrangian of fixed constraints
	for(size_t j = 0; j < n_fix_con_; j++)
	{	size_t ell        = 2 * fix_likelihood_nabs_ + j;
		w_fix_con_tmp_[j] = lambda[ell];
# if CPPAD_MIXED_HIDE_IPOPT_SCALING
		w_fix_con_tmp_[j] = scale_g_[ell] * lambda[ell];
# endif
	}
	s_vector fix_con_hes_row = mixed_object_.fix_con_hes_.subset.row();
	s_vector fix_con_hes_col = mixed_object_.fix_con_hes_.subset.col();
	d_vector fix_con_hes_val = mixed_object_.fix_con_hes_.subset.val();
	mixed_object_.fix_con_hes(
		fixed_tmp_,
		w_fix_con_tmp_,
		fix_con_hes_row,
		fix_con_hes_col,
		fix_con_hes_val
	);
	for(size_t k = 0; k < fix_con_hes_row.size(); k++)
	{	size_t index = fix_con_hes_2_lag_[k];
		assert( index < nnz_h_lag_ );
		values[index] += Number( fix_con_hes_val[k] );
	}
	//
	for(size_t ell = 0; ell < nnz_h_lag_; ell++)
	{	if( CppAD::isnan( values[ell] ) ) throw CppAD::mixed::exception(
			"", "Hessian of Lagragian has a nan"
		);
	}
	assert( size_t(nele_hess) == nnz_h_lag_ );
	return;
}
} } // END_CPPAD_MIXED_NAMESPACE
