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

$section Adaptive Step Size check of Derivative Calculations$$

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
} } // END_CPPAD_MIXED_NAMESPACE
