// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/ipopt_fixed.hpp>
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/configure.hpp>


namespace {

	// --------------------------------------------------------------------
	bool check_in_limits(double lower, double x, double upper, double tol)
	{	bool flag = true;
		if( upper >= 0.0 )
			flag &= x <= (1.0 + tol) * upper;
		else
			flag &= x <= (1.0 - tol) * upper;
		//
		if( lower >= 0.0 )
			flag &= (1.0 - tol) * lower <= x;
		else
			flag &= (1.0 + tol) * lower <= x;
		//
		return flag;
	}
}

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE


/* $$
------------------------------------------------------------------------------
$begin ipopt_fixed_eval_g$$
$spell
	CppAD
	ran_obj
	cppad
	obj
	ipopt
	bool
	const
	eval
$$

$section Compute Value of Constraint Functions$$

$head Syntax$$
$icode%ok% = eval_g(%n%, %x%, %new_x%, %m%, %g%)%$$

$head n$$
is the number of variables in the problem (dimension of x).

$head x$$
is the value for the primal variables at which the constraints
$latex g(x)$$ is computed (has size $icode n$$).

$head new_x$$
if true, no Ipopt evaluation method was previous called with the same
value for $icode x$$.

$head m$$
is the number of constraints in the problem (dimension of g(x)).

$head g$$
is set to the value for the constraint functions (has size $icode m$$).

$head ok$$
if set to false, the optimization will terminate with status set to
$cref/USER_REQUESTED_STOP
	/ipopt_fixed_finalize_solution/status/USER_REQUESTED_STOP/$$.

$head Prototype$$
$srccode%cpp% */
bool ipopt_fixed::eval_g(
	Index           n        ,  // in
	const Number*   x        ,  // in
	bool            new_x    ,  // in
	Index           m        ,  // in
	Number*         g        )  // out
/* %$$
$end
*/
{	try
	{	try_eval_g(n, x, new_x, m, g);
	}
	catch(const CppAD::mixed::exception& e)
	{	error_message_ = e.message("ipopt_fixed::eval_g");
		for(size_t j = 0; j < n_fixed_; j++)
			error_fixed_[j] = x[j];
		return false;
	}
	return true;
}
void ipopt_fixed::try_eval_g(
	Index           n        ,  // in
	const Number*   x        ,  // in
	bool            new_x    ,  // in
	Index           m        ,  // in
	Number*         g        )  // out
{
	assert( n > 0 );
	assert( size_t(n) == n_fixed_ + fix_likelihood_nabs_ );
	assert( m >= 0 );
	assert( size_t(m) == 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_ );
	//
	// fixed effects
	for(size_t j = 0; j < n_fixed_; j++)
		fixed_tmp_[j] = double( x[j] );
	//
	// check if this is a new x
	if( new_x && n_random_ > 0 )
		new_random(fixed_tmp_);
	//
	// fixed likelihood
	// (2DO: cache fix_likelihood_vec_tmp_ for eval_f with same x)
	fix_likelihood_vec_tmp_ = mixed_object_.fix_like_eval(fixed_tmp_);
	//
	// constraint part of absolute value terms in fixed likelihood
	for(size_t j = 0; j < fix_likelihood_nabs_; j++)
	{	// x[n_fixed_ + j] >= fix_likelihood_vec_tmp_[1 + j];
		assert( 2 * j + 1 < size_t(m) );
		// g[2 * j] >= 0
		g[2 * j] = Number(x[n_fixed_ + j] - fix_likelihood_vec_tmp_[1 + j]);
		// x[n_fixed_ + j] >= - fix_likelihood_vec_tmp_[1 + j]
		// g[2 * j + 1] >= 0
		g[2*j+1] = Number(x[n_fixed_ + j] + fix_likelihood_vec_tmp_[1 + j]);
	}
	//
	// fixed constraints
	assert( c_vec_tmp_.size() == n_fix_con_ );
	c_vec_tmp_ = mixed_object_.fix_con_eval(fixed_tmp_);
	for(size_t j = 0; j < n_fix_con_; j++)
	{	assert( 2 * fix_likelihood_nabs_ + j < size_t(m) );
		g[2 * fix_likelihood_nabs_ + j] = c_vec_tmp_[j];
	}
	//
	// random constraints
	assert( A_uhat_tmp_.size() == n_ran_con_ );
	mixed_object_.ran_con_eval(random_cur_, A_uhat_tmp_);
	for(size_t j = 0; j < n_ran_con_; j++)
	{	assert( 2 * fix_likelihood_nabs_ + n_fix_con_ + j < size_t(m) );
		g[2 * fix_likelihood_nabs_ + n_fix_con_ + j] = A_uhat_tmp_[j];
	}
	for(size_t i = 0; i < size_t(m); i++)
	{	if( CppAD::isnan( g[i] ) ) throw CppAD::mixed::exception(
			"", "constaint function has a nan"
		);
# if CPPAD_MIXED_HIDE_IPOPT_SCALING
		g[i] *= scale_g_[i];
# endif
	}
	return;
}
/*
-------------------------------------------------------------------------------
$begin ipopt_fixed_eval_jac_g$$
$spell
	CppAD
	ran_obj
	cppad
	obj
	ipopt
	bool
	eval
	const
	nele_jac
	Jacobian
	nnz
$$

$section Compute Jacobian of Constraint Functions$$

$head Syntax$$
$icode%ok% = eval_jac_g(
	%n%, %x%, %new_x%, %m%, %nele_jac%, %iRow%, %jCol%, %values%
)%$$

$head n$$
is the number of variables in the problem (dimension of x).

$head x$$
is the value for the primal variables at which the Jacobian
of the constraints $latex \nabla g(x)$$ is computed (has size $icode n$$).

$head new_x$$
if true, no Ipopt evaluation method was previous called with the same
value for $icode x$$.

$head m$$
is the number of constraints in the problem (dimension of g(x)).

$head nele_jac$$
is the number of non-zero elements in the Jacobian of $icode g(x)$$; i.e.,
the same as
$cref/nnz_jac_g/ipopt_fixed_get_nlp_info/nnz_jac_g/$$.

$head iRow$$
If $icode values$$ is $code NULL$$,
$icode iRow$$ has size $icode nele_jac$$ and is set to the
row indices for the non-zero entries in the Jacobian of the constraints
$latex g_x (x)$$.

$head jCol$$
If $icode values$$ is $code NULL$$,
$icode jCol$$ has size $icode nele_jac$$ and is set to the
column indices for the non-zero entries in the Jacobian of the constraints
$latex g_x (x)$$.

$head values$$
If $icode values$$ is not $code NULL$$,
it has size $icode nele_jac$$ and $icode%values%[%k%]%$$
is set to the value of element of the Jacobian $latex g_x (x)$$
with row index $icode%iRow%[%k%]%$$
and column index $icode%jCol%[%k%]%$$.

$head ok$$
if set to false, the optimization will terminate with status set to
$cref/USER_REQUESTED_STOP
	/ipopt_fixed_finalize_solution/status/USER_REQUESTED_STOP/$$.

$head Prototype$$
$srccode%cpp% */
bool ipopt_fixed::eval_jac_g(
	Index           n        ,  // in
	const Number*   x        ,  // in
	bool            new_x    ,  // in
	Index           m        ,  // in
	Index           nele_jac ,  // in
	Index*          iRow     ,  // out
	Index*          jCol     ,  // out
	Number*         values   )  // out
/* %$$
$end
*/
{	try
	{	try_eval_jac_g(n, x, new_x, m, nele_jac, iRow, jCol, values);
	}
	catch(const CppAD::mixed::exception& e)
	{	error_message_ = e.message("ipopt_fixed::eval_jac_g");
		for(size_t j = 0; j < n_fixed_; j++)
			error_fixed_[j] = x[j];
		return false;
	}
	return true;
}
void ipopt_fixed::try_eval_jac_g(
	Index           n        ,  // in
	const Number*   x        ,  // in
	bool            new_x    ,  // in
	Index           m        ,  // in
	Index           nele_jac ,  // in
	Index*          iRow     ,  // out
	Index*          jCol     ,  // out
	Number*         values   )  // out
{
	assert( n > 0 );
	assert( size_t(n) == n_fixed_ + fix_likelihood_nabs_ );
	assert( m >= 0 );
	assert( size_t(m) == 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_ );
	assert( size_t(nele_jac) == nnz_jac_g_ );
	//
	s_vector fix_like_jac_row = mixed_object_.fix_like_jac_.subset.row();
	s_vector fix_like_jac_col = mixed_object_.fix_like_jac_.subset.col();
	d_vector fix_like_jac_val = mixed_object_.fix_like_jac_.subset.val();
	s_vector fix_con_jac_row = mixed_object_.fix_con_jac_.subset.row();
	s_vector fix_con_jac_col = mixed_object_.fix_con_jac_.subset.col();
	d_vector fix_con_jac_val = mixed_object_.fix_con_jac_.subset.val();
	if( values == NULL )
	{	// just return row and column indices for l1 constraints
		size_t ell = 0;
		for(size_t k = 0; k < fix_like_jac_row.size(); k++)
		{	if( fix_like_jac_row[k] != 0 )
			{	assert( ell + 1 < nnz_jac_g_ );
				iRow[ell] = Index( 2 * fix_like_jac_row[k] - 2 );
				jCol[ell] = Index( fix_like_jac_col[k] );
				ell++;
				iRow[ell] = Index( 2 * fix_like_jac_row[k] - 1 );
				jCol[ell] = Index( fix_like_jac_col[k] );
				ell++;
			}
		}
		// auxillary variables for l1 constraints
		for(size_t j = 0; j < fix_likelihood_nabs_; j++)
		{	assert( ell + 1 < nnz_jac_g_ );
			iRow[ell] = Index( 2 * j );
			jCol[ell] = Index( n_fixed_ + j);
			ell++;
			iRow[ell] = Index( 2 * j + 1);
			jCol[ell] = Index(n_fixed_ + j);
			ell++;
		}
		// fixed constraints
		size_t offset = 2 * fix_likelihood_nabs_;
		for(size_t k = 0; k < fix_con_jac_row.size(); k++)
		{	assert( ell < nnz_jac_g_ );
			iRow[ell] = Index( offset + fix_con_jac_row[k] );
			jCol[ell] = Index( fix_con_jac_col[k] );
			ell++;
		}
		// random constraints
		offset = 2 * fix_likelihood_nabs_ + n_fix_con_;
		for(size_t k = 0; k < ran_con_jac_rcv_.nnz(); k++)
		{	assert( ell < nnz_jac_g_ );
			iRow[ell] = Index( offset + ran_con_jac_rcv_.row()[k] );
			jCol[ell] = Index( ran_con_jac_rcv_.col()[k] );
			ell++;
		}
		assert( ell == nnz_jac_g_ );
		//
		return;
	}
	assert( jac_g_row_.size() == size_t(nnz_jac_g_) );
	//
	// fixed effects
	for(size_t j = 0; j < n_fixed_; j++)
		fixed_tmp_[j] = double( x[j] );
	//
	// check if this is a new x
	if( new_x && n_random_ > 0 )
		new_random(fixed_tmp_);
	//
	// Jacobian of fixed effects likelihood
	// (2DO: do not revaluate when eval_grad_f had same x)
	mixed_object_.fix_like_jac(
		fixed_tmp_,
		fix_like_jac_row,
		fix_like_jac_col,
		fix_like_jac_val
	);
	size_t ell = 0;
	for(size_t k = 0; k < fix_like_jac_row.size(); k++)
	{	if( fix_like_jac_row[k] != 0 )
		{	assert( ell + 1 < nnz_jac_g_ );
			values[ell] = Number( - fix_like_jac_val[k] );
			ell++;
			values[ell] = Number( + fix_like_jac_val[k] );
			ell++;
		}
	}
	for(size_t j = 0; j < fix_likelihood_nabs_; j++)
	{	assert( ell + 1 < nnz_jac_g_ );
		values[ell+1] = values[ell] = Number(1.0);
		ell += 2;
	}
	//
	// Jacobian of fixed constraints
	mixed_object_.fix_con_jac(
		fixed_tmp_,
		fix_con_jac_row,
		fix_con_jac_col,
		fix_con_jac_val
	);
	for(size_t k = 0; k < fix_con_jac_row.size(); k++)
	{	assert( ell < nnz_jac_g_ );
		values[ell++] = Number( fix_con_jac_val[k] );
	}
	//
	// Jacobian of random constraints
	if( n_ran_con_ > 0 )
	{	assert( n_random_ > 0 );
		mixed_object_.ran_con_jac(fixed_tmp_, random_cur_, ran_con_jac_rcv_);
		for(size_t k = 0; k < ran_con_jac_rcv_.nnz(); k++)
		{	assert( ell < nnz_jac_g_ );
			values[ell++] = Number( ran_con_jac_rcv_.val()[k] );
		}
	}
	assert( ell == nnz_jac_g_ );
	//
	for(ell = 0; ell < nnz_jac_g_; ell++)
	{	if( CppAD::isnan( values[ell] ) ) throw CppAD::mixed::exception(
			"", "constraint Jacobian has a nan"
		);
# if CPPAD_MIXED_HIDE_IPOPT_SCALING
		size_t i     = jac_g_row_[ell];
		values[ell] *= scale_g_[i];
# endif
	}
	return;
}
/*
-------------------------------------------------------------------------------
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
/*
-------------------------------------------------------------------------------
$begin ipopt_fixed_finalize_solution$$
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
	ip
	cq
	namespace
	infeasibility
	doesn't
	Inf
	naninf
	std
	cout
	endl
	tol
	solution solution
$$

$section Get Solution Results$$

$head Syntax$$
$codei%finalize_solution(
	%status%, %n%, %x%, %z_L%, %z_U%, %m%, %g%,%$$
$icode%lambda%, %obj_value%, %ip_data%, %ip_cq%
)%$$

$head solution_$$
This routine checks the solution values and sets the member variable
$codei%
	CppAD::mixed fixed_solution solution_
%$$.

$head n$$
is the number of variables in the problem (dimension of x).

$head x$$
is the final value (best value found) for the primal variables
(has size $icode n$$).

$head z_L$$
is the final value for the $icode x$$ lower bound constraint multipliers
(has size $icode n$$).

$head z_U$$
is the final value for the $icode x$$ upper bound constraint multipliers
(has size $icode n$$).

$head m$$
is the number of constraints in the problem (dimension of g(x)).

$head lambda$$
is the final value for the g(x) constraint multipliers $latex \lambda$$.

$head obj_value$$
is the value of the objective f(x) at the final $icode x$$ value.

$head ip_data$$
Unspecified; i.e., not part of the Ipopt user API.

$head ip_cq$$
Unspecified; i.e., not part of the Ipopt user API.

$head status$$
These status values are in the $code Ipopt$$ namespace; e.g.,
$code SUCCESS$$ is short for $code Ipopt::SUCCESS$$:

$subhead SUCCESS$$
Algorithm terminated successfully at a locally optimal point,
satisfying the convergence tolerances (can be specified by options).

$subhead MAXITER_EXCEEDED$$
Maximum number of iterations exceeded (can be specified by an option).

$subhead CPUTIME_EXCEEDED$$
Maximum number of CPU seconds exceeded (can be specified by an option).

$subhead STOP_AT_TINY_STEP$$
Algorithm proceeds with very little progress.

$subhead STOP_AT_ACCEPTABLE_POINT$$
Algorithm stopped at a point that was converged, not to desired
tolerances, but to acceptable tolerances (see the acceptable-... options).

$subhead LOCAL_INFEASIBILITY$$
Algorithm converged to a point of local infeasibility. Problem may be
infeasible.

$subhead USER_REQUESTED_STOP$$
A user call-back function returned false, i.e.,
the user code requested a premature termination of the optimization.

$subhead DIVERGING_ITERATES$$
It seems that the iterates diverge.

$subhead RESTORATION_FAILURE$$
Restoration phase failed, algorithm doesn't know how to proceed.

$subhead ERROR_IN_STEP_COMPUTATION$$
An unrecoverable error occurred while Ipopt tried to compute
the search direction.

$subhead INVALID_NUMBER_DETECTED$$
Algorithm received an invalid number (such as NaN or Inf) from
the NLP; see also option check_derivatives_for_naninf.

$head Prototype$$
$srccode%cpp% */
void ipopt_fixed::finalize_solution(
	Ipopt::SolverReturn               status    ,  // in
	Index                             n         ,  // in
	const Number*                     x         ,  // in
	const Number*                     z_L       ,  // in
	const Number*                     z_U       ,  // in
	Index                             m         ,  // in
	const Number*                     g         ,  // in
	const Number*                     lambda    ,  // in
	Number                            obj_value ,  // in
	const Ipopt::IpoptData*           ip_data   ,  // in
	Ipopt::IpoptCalculatedQuantities* ip_cq     )  // in
/* %$$
$end
*/
{	bool ok = true;
	//
	assert( n > 0 );
	assert( size_t(n) == n_fixed_ + fix_likelihood_nabs_ );
	assert( m >= 0 );
	assert( size_t(m) == 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_ );
	//
	// solution_.fixed_opt
	assert( solution_.fixed_opt.size() == 0 );
	solution_.fixed_opt.resize(n_fixed_);
	for(size_t j = 0; j < n_fixed_; j++)
		solution_.fixed_opt[j] = x[j];
	//
	// solution_.fixed_lag (see below)
	//
	// solution_.ran_con_lag
	assert( solution_.ran_con_lag.size() == 0 );
	solution_.ran_con_lag.resize(n_ran_con_);
	size_t offset = 2 * fix_likelihood_nabs_ + n_fix_con_;
	for(size_t j = 0; j < n_ran_con_; j++)
		solution_.ran_con_lag[j] = lambda[ offset + j];
	//
	// short name for fixed effects tolerance
	double tol = fixed_tolerance_;
	//
	// check that x is within its limits
	for(size_t j = 0; j < n_fixed_; j++)
	{	ok &= check_in_limits(
			fixed_lower_[j], x[j], fixed_upper_[j], 2.0 * tol
		);
	}
	//
	// check that the bound multipliers are feasible
	for(size_t j = 0; j < n_fixed_ + fix_likelihood_nabs_; j++)
	{	ok &= 0.0 <= z_L[j];
		ok &= 0.0 <= z_U[j];
	}
	//
	// fixed_opt is an alias for solution_.fixed_opt
	d_vector& fixed_opt = solution_.fixed_opt;
	//
	// fixed likelihood at the final fixed effects vector
	if( fix_likelihood_vec_tmp_.size() == 0 )
		assert( mixed_object_.fix_like_eval(fixed_opt).size() == 0 );
	else
	{	fix_likelihood_vec_tmp_ = mixed_object_.fix_like_eval(fixed_opt);
		assert( fix_likelihood_vec_tmp_.size() == 1 + fix_likelihood_nabs_ );

		// check constraints corresponding to l1 terms
		for(size_t j = 0; j < fix_likelihood_nabs_; j++)
		{	double var  = double( x[n_fixed_ + j] );
			double diff = var - fix_likelihood_vec_tmp_[j + 1];
			ok         &= scale_g_[2 * j] * diff + 1e2 * tol >= 0;
			diff        = var + fix_likelihood_vec_tmp_[j + 1];
			ok         &= scale_g_[2 * j] * diff + 1e2 * tol >= 0;
		}
	}
	//
	// explicit constraints at the final fixed effects vector
	c_vec_tmp_ = mixed_object_.fix_con_eval(fixed_opt);
	assert( c_vec_tmp_.size() == n_fix_con_ );

	// check explicit constraints and set solution_.fix_con_lag
	assert( solution_.fix_con_lag.size() == 0 );
	solution_.fix_con_lag.resize(n_fix_con_);
	offset     = 2 * fix_likelihood_nabs_;
	double inf = std::numeric_limits<double>::infinity();
	for(size_t j = 0; j < n_fix_con_; j++)
	{	// It seems from testing that Ipopt is insuring the constraint to
		// be within tol, but it seems to Brad is should be within
		// scale_g_[offset +j] * tol.
		ok &= check_in_limits(
			fix_constraint_lower_[j], c_vec_tmp_[j], fix_constraint_upper_[j],
			2.0 * tol
		);
		double lam_j  = lambda[offset + j];
		double scale = 0.0;;
		if( fix_constraint_lower_[j] != -inf )
			scale = std::fabs( fix_constraint_lower_[j] );
		if( fix_constraint_upper_[j] != inf )
			scale = std::max(scale, std::fabs( fix_constraint_upper_[j] ) );
		if( scale == 0.0 )
		{	// both limits are infinity
			lam_j = 0.0;
		}
		if( c_vec_tmp_[j] - fix_constraint_lower_[j] < tol * 10. * scale )
			lam_j = std::min(lam_j, 0.0);
		if( fix_constraint_upper_[j] - c_vec_tmp_[j] < tol * 10. * scale )
			lam_j = std::max(lam_j, 0.0);
		//
		solution_.fix_con_lag[j] = lam_j;
	}
	// Evaluate gradient of f w.r.t x
	CppAD::vector<Number> grad_f(n);
	bool new_x = true;
	eval_grad_f(n, x, new_x, grad_f.data() );

	// Evaluate gradient of g w.r.t x
	CppAD::vector<Index> iRow(nnz_jac_g_), jCol(nnz_jac_g_);
	iRow.data();
	eval_jac_g(
		n, x, new_x, m, Index(nnz_jac_g_),
		iRow.data(), jCol.data(), NULL
	);
	CppAD::vector<Number> jac_g(nnz_jac_g_);
	eval_jac_g(
		n, x, new_x, m, Index(nnz_jac_g_),
		iRow.data(), jCol.data(), jac_g.data()
	);

	// Check the partial of the Lagrangian w.r.t fixed effects
	// and set solution_.fixed_lag
	assert( solution_.fixed_lag.size() == 0 );
	solution_.fixed_lag.resize(n_fixed_);
	double average = 0.0;
	for(size_t j = 0; j < n_fixed_; j++)
	{	Number sum = grad_f[j];
		for(size_t k = 0; k < nnz_jac_g_; k++)
		{	if( jCol[k] == Index(j) )
			{	Index  i = iRow[k];
				sum   += lambda[i] * jac_g[k];
			}
		}
		// sum += z_U[j] - z_L[j]; does not work because
		// Ipopt does not seem to set z_U[j] and z_L[j] accuractely
		//
		// initialize
		solution_.fixed_lag[j] = 0.0;
		//
		// scale
		double scale = std::fabs( x[j] );
		if( fixed_lower_[j] != - inf )
			scale = std::max(scale, std::fabs( fixed_lower_[j] ) );
		if( fixed_upper_[j] != + inf )
			scale = std::max(scale, std::fabs( fixed_upper_[j] ) );
		//
		// at_lower
		bool at_lower = x[j] - fixed_lower_[j] <= scale * 10. * tol;
		// catch special case where lower and upper limits are zero
		at_lower     |= fixed_lower_[j] == fixed_upper_[j];
		at_lower     &= sum > 0.0;
		if( at_lower )
			solution_.fixed_lag[j] = - sum;
		//
		// at_upper
		bool at_upper = fixed_upper_[j] - x[j] <= scale * 10. * tol;
		at_upper     |= fixed_lower_[j] == fixed_upper_[j];
		at_upper     &= sum < 0.0;
		if( at_upper )
			solution_.fixed_lag[j] = - sum;
		//
		if( ! (at_lower || at_upper) )
		{	double check;
			check = std::fabs(sum);
			if( sum >= 0.0  && fixed_lower_[j] > -inf )
				check = sum * (x[j] - fixed_lower_[j]);
			else if( sum <= 0.0 && fixed_upper_[j] < +inf )
				check = - sum * (fixed_upper_[j] - x[j]);
			assert( check >= 0.0 );
			//
			average += check / double(n_fixed_);
		}
	}
	ok &= average * scale_f_ <= 10. * tol;

	// Check the partial of the Lagrangian w.r.t auxillary variables
	average = 0.0;
	for(size_t j = n_fixed_; j < n_fixed_ + fix_likelihood_nabs_; j++)
	{	Number sum = grad_f[j];
		for(size_t k = 0; k < nnz_jac_g_; k++)
		{	if( jCol[k] == Index(j) )
			{	Index  i = iRow[k];
				sum   += lambda[i] * jac_g[k];
			}
		}
		average += std::fabs(sum) / double(n_fixed_);
	}
	ok &= average * scale_f_ <= tol;

	// set member variable finalize_solution_ok_
	finalize_solution_ok_ = ok;
}
// ---------------------------------------------------------------------------
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
/*
-------------------------------------------------------------------------------
$begin ipopt_fixed_adaptive_derivative_check$$
$spell
	CppAD
	cppad
	Ipopt
	eval
	tol
	bool
	jac
	std
$$

$section Adaptive Step Size check of eval_grad_f and eval_jac_g$$

$head Syntax$$
$icode%ok% = adaptive_derivative_check(%trace%, %relative_tol%)%$$

$head trace$$
If true, a trace of this computation is printed on standard output.

$head relative_step$$
For an unspecified set of relative step sizes between
$code 1e-3$$ and $code 1e-10$$:
If the upper and lower bounds are finite,
the step is relative to the upper minus the lower bound
(for each component of $icode x$$).
If the upper or lower bound is infinite,
the step is the maximum of the relative step
and the relative step times the absolute component of $icode x$$.

$head relative_tol$$
This is the relative tolerance for the difference between a finite difference
approximation and the evaluated derivative.
The absolute tolerance is the relative tolerance times the
sum of sum of the absolute value of the gradient and the approximation.
In the case of the Hessian, for each column the absolute value of
the diagonal element for that column is added to the sum before multiplying
by the relative tolerance.

$subhead infinity$$
If $icode trace$$ is false and $icode relative_tol$$ is equal to
$code std::numeric_limits<double>::infinity()$$,
the finite difference approximations for the derivatives are not calculated.

$head ok$$
If it is true, no function evaluation error occurred and
all the differences between the finite difference approximation
and the evaluated derivative are within the relative tolerance
(for one of the relative steps).
Otherwise,
at least one of the differences is not within the specified tolerance.

$head error_message_$$
Use $code clear_error_message()$$ to set this to the empty
string before calling $code adaptive_derivative_check$$.
(Note that it is empty after the $code ipopt_fixed$$ constructor is called.)
If upon return, $code get_error_message()$$ is non-empty,
a description of a function evaluation
error that occurred during this routine.

$head adaptive_called_$$
This member variable has prototype
$codei%
	bool adaptive_called_
%$$
It's value upon call must be $code false$$.
It is set to true at the beginning of $code adaptive_derivative_check$$,
before any other $code ipopt_fixed$$ routine is called.

$head scale_f_$$
This member variable has prototype
$codei%
	double scale_f_
%$$
and its input value does not matter.
If $icode ok$$ is true, upon return $icode scale_f_$$
is a scale factor of $latex f(x)$$; i.e., multiplier for f.
Components of $latex x$$ for which the lower and upper limits are equal
are not included in the scaling of $latex f(x)$$.

$head scale_g_$$
This member variable has prototype
$codei%
	d_vector scale_g_
%$$
and its input size is zero.
If $icode ok$$ is true, upon return $icode scale_g_$$
has size equal to range for $latex g(x)$$; i.e., $icode m$$.
For each $latex i$$, $icode%scale_g_%[%i%]%$$ is scale factor for
$latex g_i (x)$$.
Components of $latex x$$ for which the lower and upper limits are equal
are not included in the scaling of $latex f(x)$$.

$head jac_g_row_, jac_g_col_$$
These member variables have prototype
$codei%
	s_vector jac_g_row_, jac_g_col_
%$$
This input size for these vectors must be zero.
Upon return these vectors
map return index for $code eval_jac_g$$ values vector
to row and column index in $latex g'(x)$$.


$head Prototype$$
$srccode%cpp% */
bool ipopt_fixed::adaptive_derivative_check(bool trace, double relative_tol)
/* %$$
$end
*/
{	assert( error_message_ == "" );
	assert( adaptive_called_ == false );
	adaptive_called_ = true;
	//
	using std::fabs;
	// ---------------------------------------------------------------------
	// some constants
	// ---------------------------------------------------------------------
	// number of components in x
	const size_t n  = n_fixed_ + fix_likelihood_nabs_;
	// number of components if g
	const size_t m  = 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_;
	// infinity
	const double infinity  = std::numeric_limits<double>::infinity();
	// maximum scaling factor
	const double scale_max = 1e+14;
	// minimum scaling factor
	double scale_min       = 1.0 / scale_max;
	// log of maximum relatives steps size in finite differences
	double log_max_rel_step = std::log(1e-3);
	// log of minimum relative steps size in finite differences
	double log_min_rel_step = std::log(1e-10);
	// number of finite difference steps to try
	size_t n_try = 5;
	// -----------------------------------------------------------------------
	// Set scale_f_, scale_g_, to identity mapping during this routine
	// and set to to its final value just before returning.
	assert( scale_g_.size() == 0 );
	scale_f_ = 1.0;
	scale_g_.resize(m);
	for(size_t i = 0; i < m; i++)
		scale_g_[i] = 1.0;
	// -----------------------------------------------------------------------
	// Set x = x_scale
	d_vector x_scale(n);
	Number*  x = x_scale.data();
	//
	// fixed likelihood at the fixed effects scaling vector
	if( fix_likelihood_vec_tmp_.size() == 0 )
		assert( mixed_object_.fix_like_eval(fixed_scale_).size() == 0 );
	else
	{	fix_likelihood_vec_tmp_ = mixed_object_.fix_like_eval(fixed_scale_);
		assert( fix_likelihood_vec_tmp_.size() == 1 + fix_likelihood_nabs_ );
	}
	// use scaling values for fixed effects
	for(size_t j = 0; j < n_fixed_; j++)
		x_scale[j] = fixed_scale_[j];
	//
	// set auxillary variables to corresponding minimum feasible value
	for(size_t j = 0; j < fix_likelihood_nabs_; j++)
		x_scale[n_fixed_ + j] = std::fabs( fix_likelihood_vec_tmp_[1 + j] );
	// ------------------------------------------------------------------------
	// Set new_x, grad_f
	d_vector grad_f(n);
	{	bool new_x    = true;
		bool ok       = eval_grad_f( Index(n), x, new_x, grad_f.data() );
		if( ! ok )
		{	assert( error_message_ != "" );
			return false;
		}
	}
	// ------------------------------------------------------------------------
	// Set jac_g, jac_g_row_, jac_g_col_
	d_vector jac_g(nnz_jac_g_);
	{
		// get iRow and jCol for eval_jac_g
		CppAD::vector<Index> iRow(nnz_jac_g_), jCol(nnz_jac_g_);
		bool new_x     = false;
		Index nele_jac = Index( nnz_jac_g_ );
		bool ok = eval_jac_g(
			Index(n),
			x,
			new_x,
			Index(m),
			nele_jac,
			iRow.data(),
			jCol.data(),
			NULL
		);
		if( ! ok )
		{	assert( error_message_ != "" );
			return false;
		}
		assert( jac_g_row_.size() == 0 && jac_g_col_.size() == 0 );
		jac_g_row_.resize(nnz_jac_g_);
		jac_g_col_.resize(nnz_jac_g_);
		for(size_t k = 0; k < nnz_jac_g_; k++)
		{	jac_g_row_[k] = size_t( iRow[k] );
			jac_g_col_[k] = size_t( jCol[k] );
		}
		// eval_jac_g
		new_x           = false;
		nele_jac        = Index( nnz_jac_g_ );
		Number* values  = jac_g.data();
		ok = eval_jac_g(
			Index(n),
			x,
			new_x,
			Index(m),
			nele_jac,
			NULL,
			NULL,
			values
		);
		if( ! ok )
		{	assert( error_message_ != "" );
			return false;
		}
	}
	// ------------------------------------------------------------------------
	// set x_lower, x_upper, g_lower, g_upper
	d_vector x_lower(n), x_upper(n), g_lower(m), g_upper(m);
	bool ok = get_bounds_info(
		Index(n),
		x_lower.data(),
		x_upper.data(),
		Index(m),
		g_lower.data(),
		g_upper.data()
	);
	assert( ok );
	// ------------------------------------------------------------------------
	// max_jac_g : maximum absolute partial for each component of g
	d_vector max_jac_g(m);
	for(size_t i = 0; i < m; i++)
		max_jac_g[i] = scale_min;
	for(size_t k = 0; k < nnz_jac_g_; k++)
	{	size_t i = jac_g_row_[k];
		size_t j = jac_g_col_[k];
		if( x_lower[j] < x_upper[j] )
			max_jac_g[i] = std::max( max_jac_g[i], std::fabs( jac_g[k] ) );
	}
	// -----------------------------------------------------------------------
	// max_grad_f : maximum absolute partial for objective.
	// Skip absolute value partials in f, which are one, but
	// include the corresponding partials in g
	double max_grad_f = scale_min;
	//
	// fixed effect terms
	for(size_t j = 0; j < n_fixed_; j++) if( x_lower[j] < x_upper[j] )
		max_grad_f = std::max( max_grad_f, std::fabs( grad_f[j] ) );
	//
	// skip axuillary variable terms (gradients are one)
# ifndef NDEBUG
	for(size_t j = 0; j < fix_likelihood_nabs_; j++)
		assert( grad_f[n_fixed_ + j] == 1.0 );
# endif
	// include absolute value terms in objective scaling
	for(size_t i = 0; i < fix_likelihood_nabs_; i++)
	{	assert( max_jac_g[2*i] == max_jac_g[2*i+1] );
		max_grad_f = std::max( max_grad_f, std::fabs( max_jac_g[2*i] ) );
	}
	// -----------------------------------------------------------------------
	// scale_f
	double scale_f = std::max( scale_min, 1.0 / max_grad_f);
	scale_f        = std::min( scale_max, scale_f);
	// ----------------------------------------------------------------------
	// scale_g
	d_vector scale_g(m);
	for(size_t i = 0; i < m; i++)
	{	if( i < 2 * fix_likelihood_nabs_ )
			scale_g[i] = scale_f;
		else
		{	scale_g[i] = std::max( scale_min, 1.0 / max_jac_g[i] );
			scale_g[i] = std::min( scale_max, scale_g[i] );
		}
	}
	// ----------------------------------------------------------------------
	// If no need to compute finite differences, set scaling and return
	if( (! trace) && relative_tol == infinity )
	{
		scale_f_       = scale_f;
		scale_g_       = scale_g;
		return true;
	}
	// ------------------------------------------------------------------------
	// Test the derivatvie
	// ------------------------------------------------------------------------
	// set obj_value = f(x)
	double obj_value;
	{
		bool new_x = false;
		ok = eval_f(Index(n), x, new_x, obj_value);
		if( ! ok )
		{	assert( error_message_ != "" );
			return false;
		}
	}
	//
	// set con_value = g(x)
	d_vector con_value(m);
	{
		bool new_x = false;
		ok = eval_g(Index(n), x, new_x, Index(m), con_value.data() );
		if( ! ok )
		{	assert( error_message_ != "" );
			return false;
		}
	}
	//
	// set hes_value = Hessian of obj_factor * f(x) + sum lambda_i * g_i (x)
	d_vector hes_value( nnz_h_lag_ );
	Number   obj_factor = 2.0;
	d_vector lambda(m);
	for(size_t i = 0; i < m; i++)
		lambda[i] = double(i + 1 + m) / double(m);
	if( mixed_object_.quasi_fixed_ == false )
	{	bool new_x      = false;
		bool new_lambda = true;
		Index nele_hes  = Index(nnz_h_lag_);
		Index*  iRow    = NULL;
		Index*  jCol    = NULL;
		ok = eval_h(
			Index(n),
			x,
			new_x,
			obj_factor,
			Index(m),
			lambda.data(),
			new_lambda,
			nele_hes,
			iRow,
			jCol,
			hes_value.data()
		);
	}
	//
	// difference of log of relative step between trys
	double log_diff = (log_max_rel_step - log_min_rel_step) / double(n_try-1);
	//
	// initialize x_step = x_scale
	d_vector x_step(x_scale);
	// ------------------------------------------------------------------------
	// check grad_f
	size_t line_count = 0;
	for(size_t j = 0; j < n; j++) if( x_lower[j] < x_upper[j] )
	{	double abs_obj         = fabs(obj_value);
		double best_err        = infinity;
		double best_step       = infinity;
		double best_approx     = infinity;
		//
		// loop over relative step sizes
		size_t i_try           = 0;
		while( i_try < n_try && best_err > relative_tol )
		{	double log_next      = log_min_rel_step + log_diff * double(i_try);
			double relative_step = std::exp(log_next);
			//
			// step size
			double step = relative_step * ( x_upper[j] - x_lower[j] );
			if( x_upper[j] == nlp_upper_bound_inf_ ||
				x_lower[j] == nlp_lower_bound_inf_  )
			{	step = relative_step * fabs( x_scale[j] );
				step = std::max( step, relative_step );
			}

			// x_plus, obj_plus
			double obj_plus;
			double x_plus = std::min(x_scale[j] + step, x_upper[j]);
			x_step[j]     = x_plus;
			bool new_x    = true;
			ok = eval_f(Index(n), x_step.data(), new_x, obj_plus);
			if( ! ok )
			{	assert( error_message_ != "" );
				return false;
			}
			abs_obj = std::max(abs_obj, fabs(obj_plus) );

			// x_minus, obj_minus
			double obj_minus;
			double x_minus = std::max(x_scale[j] - step, x_lower[j]);
			x_step[j]      = x_minus;
			new_x          = true;
			eval_f(Index(n), x_step.data(), new_x, obj_minus);
			if( ! ok )
			{	assert( error_message_ != "" );
				return false;
			}
			abs_obj = std::max(abs_obj, fabs(obj_plus) );
			abs_obj = std::max(abs_obj, fabs(obj_minus) );

			// restore j-th component of x_step
			x_step[j]      = x_scale[j];

			// finite difference approximation for derivative
			double approx = (obj_plus - obj_minus) / (x_plus - x_minus);

			// relative difference
			double diff           = grad_f[j] - approx;
			double denominator    = fabs(grad_f[j]) + fabs(approx);
			if( denominator == 0.0 )
				denominator = 1.0;
			double relative_err  = fabs(diff) / denominator;

			// best
			if( relative_err < best_err )
			{	best_err    = relative_err;
				best_step   = step;
				best_approx = approx;
			}
			//
			// next try
			++i_try;
		}
		// trace
		bool trace_j = trace || best_err > relative_tol;
		if( trace_j )
		{	if( line_count % 20 == 0 )
				std::cout << std::endl
					<< std::right
					<< std::setw(4)  << "j"
					<< std::setw(11) << "step"
					<< std::setw(11) << "f"
					<< std::setw(11) << "grad"
					<< std::setw(11) << "apx"
					<< std::setw(11) << "err"
					<< std::endl;
			std::cout
				<< std::setprecision(4)
				<< std::setw(4)  << j
				<< std::setw(11) << best_step
				<< std::setw(11) << obj_value
				<< std::setw(11) << grad_f[j]
				<< std::setw(11) << best_approx
				<< std::setw(11) << best_err
				<< std::endl;
			line_count++;
		}
		//
		// ok
		ok &= ok && best_err <= relative_tol;
	}
	// ------------------------------------------------------------------------
	// check jac_g
	line_count              = 0;
	double max_best_err_all = 0.0;
	for(size_t j = 0; j < n; j++) if( x_lower[j] < x_upper[j] )
	{	d_vector best_err(m), best_step(m), best_approx(m), jac(m);
		for(size_t i = 0; i < m; i++)
		{	best_err[i]     = infinity;
			best_step[i]    = infinity;
			best_approx[i]  = infinity;
			jac[i]          = 0.0;
		}
		// value of this coluimn of the jacobian
		for(size_t k = 0; k < nnz_jac_g_; k++) if( jac_g_col_[k] == j )
		{	size_t i = jac_g_row_[k];
			jac[i]   = jac_g[k];
		}
		//
		// loop over relative step sizes
		size_t i_try        = 0;
		double max_best_err = infinity;
		while( i_try < n_try && max_best_err > relative_tol )
		{	double log_next      = log_min_rel_step + log_diff * double(i_try);
			double relative_step = std::exp(log_next);
			//
			// step size
			double step = relative_step * ( x_upper[j] - x_lower[j] );
			if( x_upper[j] == nlp_upper_bound_inf_ ||
				x_lower[j] == nlp_lower_bound_inf_  )
			{	step = relative_step * fabs( x_scale[j] );
				step = std::max( step, relative_step );
			}

			// x_plus, con_plus
			d_vector con_plus(m);
			double x_plus = std::min(x_scale[j] + step, x_upper[j]);
			x_step[j]     = x_plus;
			bool new_x    = true;
			ok = eval_g(
				Index(n), x_step.data(), new_x, Index(m), con_plus.data()
			);
			if( ! ok )
			{	assert( error_message_ != "" );
				return false;
			}
			// x_minus con_minus
			d_vector con_minus(m);
			double x_minus = std::max(x_scale[j] - step, x_lower[j]);
			x_step[j]      = x_minus;
			new_x          = true;
			ok = eval_g(
				Index(n), x_step.data(), new_x, Index(m), con_minus.data()
			);
			if( ! ok )
			{	assert( error_message_ != "" );
				return false;
			}

			// actual step size
			step = x_plus - x_minus;

			// restore j-th component of x_step
			x_step[j]  = x_scale[j];

			max_best_err = 0.0;
			for(size_t i = 0; i < m; i++)
			{
				// finite difference approximation for derivative
				double approx = (con_plus[i] - con_minus[i])/ step;

				// relative difference
				double diff           = jac[i] - approx;
				double denominator    = fabs(jac[i]) + fabs(approx);
				if( denominator == 0.0 )
					denominator = 1.0;
				double relative_err  = fabs(diff) / denominator;

				// best
				if( relative_err < best_err[i] )
				{	best_err[i]    = relative_err;
					best_step[i]   = step;
					best_approx[i] = approx;
				}
				max_best_err = std::max(max_best_err, best_err[i]);
			}
			//
			// next try
			++i_try;
		}
		max_best_err_all = std::max(max_best_err_all, max_best_err);
		for(size_t i = 0; i < m; i++)
		{	// trace
			bool trace_ij = trace || best_err[i] > relative_tol;
			if( trace_ij )
			{	if( line_count % 20 == 0 )
					std::cout << std::endl
						<< std::right
						<< std::setw(4)  << "i"
						<< std::setw(4)  << "j"
						<< std::setw(11) << "step"
						<< std::setw(11) << "f"
						<< std::setw(11) << "grad"
						<< std::setw(11) << "apx"
						<< std::setw(11) << "err"
						<< std::endl;
				std::cout
					<< std::setprecision(4)
					<< std::setw(4)  << i
					<< std::setw(4)  << j
					<< std::setw(11) << best_step[i]
					<< std::setw(11) << con_value[i]
					<< std::setw(11) << jac[i]
					<< std::setw(11) << best_approx[i]
					<< std::setw(11) << best_err[i]
					<< std::endl;
				line_count++;
			}
		}
		//
		// ok
		ok &= ok && max_best_err_all <= relative_tol;
	}
	// ------------------------------------------------------------------------
	// check hes_value
	line_count       = 0;
	max_best_err_all = 0.0;
	if( mixed_object_.quasi_fixed_ == false )
	for(size_t j2 = 0; j2 < n; j2++)
	if( x_lower[j2] < x_upper[j2] )
	{
		d_vector best_err(n), best_step(n), best_approx(n), hess(n);
		for(size_t j1 = 0; j1 < n; j1++)
		{	best_err[j1]    = infinity;
			best_step[j1]   = infinity;
			best_approx[j1] = infinity;
			hess[j1]        = 0.0;
		}
		// value of this column of the hessian
		for(size_t k = 0; k < nnz_h_lag_; k++)
		if( lag_hes_col_[k] == j2 )
		{	size_t j1 = lag_hes_row_[k];
			hess[j1]  = hes_value[k];
		}

		//
		// loop over relative step sizes
		size_t i_try        = 0;
		double max_best_err = infinity;
		while( i_try < n_try && max_best_err > relative_tol )
		{	double log_next      = log_min_rel_step + log_diff * double(i_try);
			double relative_step = std::exp(log_next);
			//
			// step size
			double step = relative_step * ( x_upper[j2] - x_lower[j2] );
			if( x_upper[j2] == nlp_upper_bound_inf_ ||
				x_lower[j2] == nlp_lower_bound_inf_  )
			{	step = relative_step * fabs( x_scale[j2] );
				step = std::max( step, relative_step );
			}

			// x_plus, grad_f_plus, jac_g_plus
			d_vector grad_f_plus(n);
			d_vector jac_g_plus(nnz_jac_g_);
			double x_plus = std::min(x_scale[j2] + step, x_upper[j2]);
			x_step[j2]    = x_plus;
			bool new_x    = true;
			Index* iRow   = NULL;
			Index* jCol   = NULL;
			ok  = eval_grad_f(
				Index(n), x_step.data(), new_x, grad_f_plus.data()
			);
			Index nele_jac = Index( nnz_jac_g_ );
			ok &= eval_jac_g(
				Index(n),
				x_step.data(),
				new_x,
				Index(m),
				nele_jac,
				iRow,
				jCol,
				jac_g_plus.data()
			);
			if( ! ok )
			{	assert( error_message_ != "" );
				return false;
			}

			// x_minus, grad_f_minus, jac_g_minus
			d_vector grad_f_minus(n);
			d_vector jac_g_minus(nnz_jac_g_);
			double x_minus = std::min(x_scale[j2] - step, x_upper[j2]);
			x_step[j2]    = x_minus;
			new_x         = true;
			ok  = eval_grad_f(
				Index(n), x_step.data(), new_x, grad_f_minus.data()
			);
			ok &= eval_jac_g(
				Index(n),
				x_step.data(),
				new_x,
				Index(m),
				nele_jac,
				iRow,
				jCol,
				jac_g_minus.data()
			);
			if( ! ok )
			{	assert( error_message_ != "" );
				return false;
			}

			// actual step size
			step = x_plus - x_minus;

			// restore j-th component of x_step
			x_step[j2] = x_scale[j2];

			// Initailize j-th column of Hessian with f contribution
			d_vector approx(n);
			for(size_t j1 = 0; j1 < n; j1++)
			{	double d2f =(grad_f_plus[j1] - grad_f_minus[j1]) / step;
				approx[j1] = obj_factor * d2f;
			}
			// add in contribution for g
			for(size_t k = 0; k < nnz_jac_g_; k++)
			{	size_t i    = jac_g_row_[k];
				size_t j1   = jac_g_col_[k];
				double d2g   = (jac_g_plus[k] - jac_g_minus[k])/ step;
				approx[j1] += lambda[i] * d2g;
			}
			//
			// absolute value of element in this column and on diagonal
			double abs_diag = fabs( hess[j2] );
			//
			// only check the lower trinagle
			max_best_err  = 0.0;
			for(size_t j1 = j2; j1 < n; j1++)
			{	// relative difference
				double diff  = hess[j1] - approx[j1];
				double den   = abs_diag + fabs(hess[j1]) + fabs(approx[j1]);
				if( den == 0.0 )
					den = 1.0;
				double relative_err  = fabs(diff) / den;

				// best
				if( relative_err < best_err[j1] )
				{	best_err[j1]    = relative_err;
					best_step[j1]   = step;
					best_approx[j1] = approx[j1];
				}
				max_best_err = std::max(max_best_err, best_err[j1]);
			}
			//
			// next try
			++i_try;
		}
		max_best_err_all = std::max(max_best_err_all, max_best_err);
		//
		// only display the lower triangle
		for(size_t j1 = j2; j1 < n; j1++)
		{	// trace
			bool trace_j1 = trace || best_err[j1] > relative_tol;
			if( trace_j1 )
			{	if( line_count % 20 == 0 )
					std::cout << std::endl
						<< std::right
						<< std::setw(4)  << "j1"
						<< std::setw(4)  << "j2"
						<< std::setw(11) << "step"
						<< std::setw(11) << "hess"
						<< std::setw(11) << "apx"
						<< std::setw(11) << "err"
						<< std::endl;
				std::cout
					<< std::setprecision(4)
					<< std::setw(4)  << j1
					<< std::setw(4)  << j2
					<< std::setw(11) << best_step[j1]
					<< std::setw(11) << hess[j1]
					<< std::setw(11) << best_approx[j1]
					<< std::setw(11) << best_err[j1]
					<< std::endl;
				line_count++;
			}
		}
		//
		// ok
		ok &= ok && max_best_err_all <= relative_tol;
	}
	// -----------------------------------------------------------------------
	// Set scaling
	scale_f_       = scale_f;
	scale_g_       = scale_g;
	return ok;
}
/*$
$begin ipopt_fixed_new_random$$
$spell
	vec
	ipopt
	const
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
// ---------------------------------------------------------------------------
} } // END_CPPAD_MIXED_NAMESPACE
