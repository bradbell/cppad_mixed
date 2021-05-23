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

// callback used by one_dim_derivative_chk
bool ipopt_fixed::one_dim_function(double x_in, d_vector& fun_out)
{	// new_x
	bool    new_x          = true;
	//
	// m: number of components in one dimensional function range space
	size_t  m              = fun_out.size();
	//
	// n: number of arguemnts to function being optimized
	size_t  n              = n_fixed_ + fix_likelihood_nabs_;
	//
	// j: component in argument space were one dimensional function defined
	size_t  j              = one_dim_function_j_;
	//
	// x_save: original value of j-th component of x (used to restore value)
	double  x_save         = one_dim_function_x_[j];
	//
	// replace the j-th compponent with requested value
	one_dim_function_x_[j] = x_in;
	//
	// x: pointer to argument for n dimensional function
	Number* x              = one_dim_function_x_.data();
	//
	// ok: did the function evaluation complete
	bool ok = false;
	//
	switch( one_dim_function_eval_ )
	{	// evaluting f(x) along j-th component of x
		case eval_f_enum:
		ok = eval_f(Index(n), x, new_x, fun_out[0]);
		break;
		//
		// evaluating g(x) along j-th component of x
		case eval_g_enum:
		ok = eval_g(
			Index(n), x, new_x, Index(m), fun_out.data()
		);
		break;
	}
	// restore one_dimensional_function_x_
	one_dim_function_x_[j] = x_save;
	//
	// return function evaluation status
	return ok;
}


/*
$begin ipopt_fixed_adaptive_derivative_check$$
$spell
	differenced
	CppAD
	cppad
	Ipopt
	eval
	tol
	bool
	jac
	std
$$

$section Adaptive Step Size check of Derivative Calculations$$

$head Syntax$$
$icode%ok% = adaptive_derivative_check(%trace%, %relative_tol%)%$$

$head trace$$
If true, a trace of this computation is printed on standard output.

$head relative_step$$
For an unspecified set of relative step sizes between
$code 1e-1$$ and $code 1e-10$$:
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
sum of sum of the absolute value of the gradient, the approximation.
In addition, the 100 times the square root of machine epsilon time
the size of the values being differenced was added to the sum.
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

$head fixed_scale_$$
This member variable has prototype
$codei%
	d_vector fixed_scale_
%$$
and is size is equal to the number of fixed effects.
It is used as the point where in the extended space $icode x$$
where the scaling and derivative testing is done.

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
	//
	// sqrt_eps: square root of machine epsilon
	double sqrt_eps = std::sqrt( std::numeric_limits<double>::epsilon() );
	//
	// n: number of components in x
	const size_t n  = n_fixed_ + fix_likelihood_nabs_;
	//
	// m: number of components if g
	const size_t m  = 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_;
	//
	// infinity
	const double infinity  = std::numeric_limits<double>::infinity();
	//
	// scale_max: maximum scaling factor
	const double scale_max = 1e+14;
	//
	// scale_min: minimum scaling factor
	double scale_min       = 1.0 / scale_max;
	//
	// log_max_rel_step: maximum log relatives steps size in finite differences
	double log_max_rel_step = std::log(0.1);
	//
	// log_min_rel_step: minimum log relative steps size in finite differences
	double log_min_rel_step = std::log(1e-10);
	//
	// n_try: number of finite difference steps to try
	size_t n_try = 5;
	//
	// scale_f_, scale_g_: are the identity mapping during this routine
	// and set to to its final value just before returning.
	assert( scale_g_.size() == 0 );
	scale_f_ = 1.0;
	scale_g_.resize(m);
	for(size_t i = 0; i < m; i++)
		scale_g_[i] = 1.0;
	//
	// x, x_scale: two view of same vector of argument values
	d_vector x_scale(n);
	Number*  x = x_scale.data();
	// set fixed effect in x_scale
	for(size_t j = 0; j < n_fixed_; j++)
		x_scale[j] = fixed_scale_[j];
	// set auxillary variables in x_scale
	for(size_t j = 0; j < fix_likelihood_nabs_; j++)
		x_scale[n_fixed_ + j] = std::fabs( fix_likelihood_vec_tmp_[1 + j] );
	//
	// fixed_likelihood_vec_tmp_: set to value offixed likelihood at the
	// fixed effects scaling vector
	if( fix_likelihood_vec_tmp_.size() == 0 )
		assert( mixed_object_.fix_like_eval(fixed_scale_).size() == 0 );
	else
	{	fix_likelihood_vec_tmp_ = mixed_object_.fix_like_eval(fixed_scale_);
		assert( fix_likelihood_vec_tmp_.size() == 1 + fix_likelihood_nabs_ );
	}
	//
	// grad_f: gradinet of f(x)
	d_vector grad_f(n);
	{	bool new_x    = true;
		bool ok       = eval_grad_f( Index(n), x, new_x, grad_f.data() );
		if( ! ok )
		{	assert( error_message_ != "" );
			return false;
		}
	}
	//
	// jac_g, jac_g_row_, jac_g_col_: sparse representation of Jacobian of g(x)
	d_vector jac_g(nnz_jac_g_);
	{	//
		// iRow, jCol: sparsity pattern for Jacobian of g(x)
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
		//
		// jac_g_row_, jac_g_col_
		assert( jac_g_row_.size() == 0 && jac_g_col_.size() == 0 );
		jac_g_row_.resize(nnz_jac_g_);
		jac_g_col_.resize(nnz_jac_g_);
		for(size_t k = 0; k < nnz_jac_g_; k++)
		{	jac_g_row_[k] = size_t( iRow[k] );
			jac_g_col_[k] = size_t( jCol[k] );
		}
		//
		// jac_g
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
	//
	// x_lower, x_upper, g_lower, g_upper
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
	// Skip absolute value partials in f, which are one, but
	// include the corresponding partials in g
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
	// scale_f
	double scale_f = std::max( scale_min, 1.0 / max_grad_f);
	scale_f        = std::min( scale_max, scale_f);
	//
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
	//
	// If not tracing or checking derivatives, set scaling and return
	if( (! trace) && relative_tol == infinity )
	{	// scale_f_, scale_g_
		scale_f_       = scale_f;
		scale_g_       = scale_g;
		return true;
	}
	// ------------------------------------------------------------------------
	// Test the derivatvie
	// ------------------------------------------------------------------------
	// obj_value: value of the objective function f(x)
	double obj_value;
	{	bool new_x = false;
		ok = eval_f(Index(n), x, new_x, obj_value);
		if( ! ok )
		{	assert( error_message_ != "" );
			return false;
		}
	}
	//
	// con_value: value of the constraint function g(x)
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
	// obj_factor:
	Number   obj_factor = 2.0;
	//
	// lambda: Lagrange multipliers`
	d_vector lambda(m);
	for(size_t i = 0; i < m; i++)
		lambda[i] = double(i + 1 + m) / double(m);
	//
	// L(x) = obj_factor * f(x) = sum_i lambda[i] * g_i(x)
	//
	// hes_value: values in sparse representation of Hessian
	d_vector hes_value( nnz_h_lag_ );
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
	// log_diff: difference of log of relative step between trys
	double log_diff = (log_max_rel_step - log_min_rel_step) / double(n_try-1);
	//
	// initialize x_step = x_scale
	d_vector x_step(x_scale);
	// ------------------------------------------------------------------------
	// check grad_f
	size_t line_count = 0;
	one_dim_function_x_    = x_scale;
	one_dim_function_eval_ = eval_f_enum;
	d_vector obj_value_vec(1);
	obj_value_vec[0] = obj_value;
	//
	d_vector grad_f_j_vec(1);
	for(size_t j = 0; j < n; j++) if( x_lower[j] < x_upper[j] )
	{	one_dim_function_j_ = j;
		double x_low = x_lower[j];
		if( x_low == nlp_lower_bound_inf_ )
			x_low = -infinity;
		double x_up = x_upper[j];
		if( x_up == nlp_upper_bound_inf_ )
			x_up = infinity;
		grad_f_j_vec[0] = grad_f[j];
		one_dim_derivative_chk_result result = one_dim_derivative_chk(
			*this,
			x_low,
			x_up,
			x_scale[j],
			obj_value_vec,
			grad_f_j_vec,
			relative_tol
		);
		if( result.rel_err[0] == infinity )
		{	assert( error_message_ != "" );
			return false;
		}
		// trace
		bool trace_j = trace || result.rel_err[0] > relative_tol;
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
				<< std::setw(11) << result.step[0]
				<< std::setw(11) << obj_value_vec[0]
				<< std::setw(11) << grad_f_j_vec[0]
				<< std::setw(11) << result.apx_dfdx[0]
				<< std::setw(11) << result.rel_err[0]
				<< std::endl;
			line_count++;
		}
		// ok
		ok &= ok && result.rel_err[0] <= relative_tol;
	}
	// ------------------------------------------------------------------------
	// check jacobian of g
	line_count             = 0;
	double rel_err_max     = 0.0;
	one_dim_function_x_    = x_scale;
	one_dim_function_eval_ = eval_g_enum;
	for(size_t j = 0; j < n; j++) if( x_lower[j] < x_upper[j] )
	{	// jac_g_j
		d_vector jac_g_j(m);
		for(size_t i = 0; i < m; i++)
			jac_g_j[i] = 0.0;
		for(size_t k = 0; k < nnz_jac_g_; k++) if( jac_g_col_[k] == j )
		{	size_t i   = jac_g_row_[k];
			jac_g_j[i] = jac_g[k];
		}
		//
		one_dim_function_j_ = j;
		double x_low = x_lower[j];
		if( x_low == nlp_lower_bound_inf_ )
			x_low = -infinity;
		double x_up = x_upper[j];
		if( x_up == nlp_upper_bound_inf_ )
			x_up = infinity;
		one_dim_derivative_chk_result result = one_dim_derivative_chk(
			*this,
			x_low,
			x_up,
			x_scale[j],
			con_value,
			jac_g_j,
			relative_tol
		);
		// rel_err_max
		rel_err_max = 0.0;
		for(size_t i = 0; i < m; ++i)
			rel_err_max = std::max(rel_err_max, result.rel_err[i]);
		if( rel_err_max == infinity )
		{	assert( error_message_ != "" );
			return false;
		}
		for(size_t i = 0; i < m; i++)
		{	// trace
			bool trace_ij = trace || result.rel_err[i] > relative_tol;
			if( trace_ij )
			{	if( line_count % 20 == 0 )
					std::cout << std::endl
						<< std::right
						<< std::setw(4)  << "i"
						<< std::setw(4)  << "j"
						<< std::setw(11) << "step"
						<< std::setw(11) << "g[i]"
						<< std::setw(11) << "jac[i,j]"
						<< std::setw(11) << "apx"
						<< std::setw(11) << "err"
						<< std::endl;
				std::cout
					<< std::setprecision(4)
					<< std::setw(4)  << i
					<< std::setw(4)  << j
					<< std::setw(11) << result.step[i]
					<< std::setw(11) << con_value[i]
					<< std::setw(11) << jac_g_j[i]
					<< std::setw(11) << result.apx_dfdx[i]
					<< std::setw(11) << result.rel_err[i]
					<< std::endl;
				line_count++;
			}
		}
		// ok
		ok &= ok && rel_err_max <= relative_tol;
	}
	// ------------------------------------------------------------------------
	// check hes_value
	line_count       = 0;
	double max_best_err_all = 0.0;
	if( mixed_object_.quasi_fixed_ == false )
	for(size_t j2 = 0; j2 < n; j2++)
	if( x_lower[j2] < x_upper[j2] )
	{
		d_vector best_err(n), best_step(n), best_approx(n), best_grad_L(n),hess(n);
		for(size_t j1 = 0; j1 < n; j1++)
		{	best_err[j1]    = infinity;
			best_step[j1]   = infinity;
			best_approx[j1] = infinity;
			best_grad_L[j1] = infinity;
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
		{	double log_next      = log_max_rel_step - log_diff * double(i_try);
			double relative_step = std::exp(log_next);
			//
			// step size
			double step = relative_step * fabs( x_scale[j2] );
			if( x_upper[j2] != nlp_upper_bound_inf_ )
				step = std::max(step, relative_step * fabs( x_upper[j2] ) );
			if( x_lower[j2] != nlp_lower_bound_inf_  )
				step = std::max(step, relative_step * fabs( x_lower[j2] ) );
			if( step == 0.0 )
				step = relative_step;

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
			d_vector grad_L(n);
			for(size_t j1 = 0; j1 < n; j1++)
			{	double d2f =(grad_f_plus[j1] - grad_f_minus[j1]) / step;
				grad_L[j1] =(grad_f_plus[j1] + grad_f_minus[j1]) / 2.0;
				approx[j1] = obj_factor * d2f;
			}
			// add in contribution for g
			for(size_t k = 0; k < nnz_jac_g_; k++)
			{	size_t i    = jac_g_row_[k];
				size_t j1   = jac_g_col_[k];
				double d2g   = (jac_g_plus[k] - jac_g_minus[k])/ step;
				approx[j1] += lambda[i] * d2g;
				grad_L[j1] += lambda[i] * (jac_g_plus[k] + jac_g_minus[k]) / 2.0;
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
				den         += 100.0 * sqrt_eps * grad_L[j1];
				double relative_err = fabs(diff);
				if( den > 0.0 )
					relative_err  = fabs(diff) / den;


				// best
				if( relative_err < best_err[j1] )
				{	best_err[j1]    = relative_err;
					best_step[j1]   = step;
					best_approx[j1] = approx[j1];
					best_grad_L[j1] = grad_L[j1];
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
						<< std::setw(11) << "grad_L"
						<< std::setw(11) << "hess"
						<< std::setw(11) << "apx"
						<< std::setw(11) << "err"
						<< std::endl;
				std::cout
					<< std::setprecision(4)
					<< std::setw(4)  << j1
					<< std::setw(4)  << j2
					<< std::setw(11) << best_step[j1]
					<< std::setw(11) << best_grad_L[j1]
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
} } // END_CPPAD_MIXED_NAMESPACE
