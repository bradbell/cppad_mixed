/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/ipopt_random.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
$begin ipopt_random_eval_h$$
$spell
	eval
	obj
	nele
	hess
	Ipopt
	nnz
	Teylor
	vec
	hes
	Taylor
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
$cref/nnz_h_lag/ipopt_random_get_nlp_info/nnz_h_lag/$$.

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
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

$head mixed_object_.ran_like_fun_$$
if $icode new_x$$ is true,
after this call, the zero order Taylor coefficients in this function
will corresponding to the value of $icode fixed_vec_$$ and
the random effects in $icode x$$.

$head mixed_object_.ran_hes_fun_$$
After this call, the zero order Taylor coefficients in this function
will corresponding to the value of $icode fixed_vec_$$ and
the random effects in $icode x$$.

$head Prototype$$
$srccode%cpp% */
bool ipopt_random::eval_h(
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
	{	error_message_ = e.message("ipopt_random::eval_h");
		for(size_t j = 0; j < n_random_; j++)
			error_random_[j] = x[j];
		return false;
	}
	return true;
}
void ipopt_random::try_eval_h(
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
{	assert( size_t(n) == n_random_ );
	assert( m == 0 );
	assert( size_t(nele_hess) == nnz_h_lag_ );
	//
	if( new_x )
	{	// set the zero order Taylor coefficients in
		// mixed_object_.ran_like_fun_
		Number obj_value;
		eval_f(n, x, new_x, obj_value);
	}
	//
	const s_vector& row( mixed_object_.ran_hes_rcv_.row() );
	const s_vector& col( mixed_object_.ran_hes_rcv_.col() );
	assert( row.size() == nnz_h_lag_ );
	assert( col.size() == nnz_h_lag_ );
	//
	if( values == NULL )
	{	for(size_t k = 0; k < nnz_h_lag_; k++)
		{	// only returning lower triagle of Hessian of objective
			assert( n_fixed_ <= col[k] );
			assert( col[k] <= row[k] );
			assert( row[k] <= n_fixed_ + n_random_ );
			iRow[k] = static_cast<Index>( row[k] - n_fixed_ );
			jCol[k] = static_cast<Index>( col[k] - n_fixed_ );
		}
		assert( ! new_x );
		return;
	}
	//
	// random effects as a vector
	d_vector random_vec(n_random_);
	for(size_t j = 0; j < n_random_; j++)
		random_vec[j] = x[j];
	//
	// pack both the fixed and random effects into one vector
	d_vector both_vec(n_fixed_ + n_random_);
	mixed_object_.pack(fixed_vec_, random_vec, both_vec);
	//
	if( new_x )
	{	// set zero order Taylor coefficient in ran_like_fun_.
		d_vector vec = mixed_object_.ran_like_fun_.Forward(0, both_vec);
		if( CppAD::hasnan( vec ) ) throw CppAD::mixed::exception(
			"", "Hessian has a nan"
		);
	}
	//
	// computes the Hessian of objecive w.r.t random effects  f_uu (theta, u)
	d_vector val = mixed_object_.ran_hes_fun_.Forward(0, both_vec);
	if( CppAD::hasnan( val ) ) throw CppAD::mixed::exception(
		"", "Hessian has a nan"
	);
	assert( val.size() == nnz_h_lag_ );
	//
	// return the values
	for(size_t k = 0; k < nnz_h_lag_; k++)
		values[k] = obj_factor * static_cast<Number>( val[k] );
	//
	return;
}
} } // END_CPPAD_MIXED_NAMESPACE
