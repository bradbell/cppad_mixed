# ifndef CPPAD_MIXED_ONE_DIM_DERIVATIVE_CHK_HPP
# define CPPAD_MIXED_ONE_DIM_DERIVATIVE_CHK_HPP

/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin one_dim_derivative_chk$$
$spell
	apx
	chk
	obj
	dfdx
	rel_tol
	bool
$$

$section One Dimensional Finite Difference Derivative Check$$

$head Syntax$$
$icode%result% = one_dim_derivative_chk(
	%obj%, %x_lower%, %x_upper%, %x%, %f%, %dfdx%, %rel_tol%
)%$$

$head Prototype$$
$srcthisfile%0%// BEGIN_PROTOTYPE%// END_PROTOTYPE%1%$$

$head obj$$
This is a object such that the syntax
$codei%
	%ok% = %obj%.one_dim_function(%x_in%, %f_out%)
%$$
evaluates the function.
The argument $icode x_in$$ has prototype
$codei%
	double x
%$$
and is the point at which we are evaluating the function.
The argument $icode f_out$$ has prototype
$codei%
	d_vector& f_out
%$$
Its input value does not matter and upon return it is the
value of the function at $icode x_in$$.
The return value $icode ok$$ has prototype
$codei%
	bool %ok%
%$$
If it is true, the computation of $icode f_out$$ completed.
Otherwise, an error was detected and the computation was aborted.

$head x_lower$$
This is a lower limit for the argument $icode x_in$$ to $icode fun$$.
If it is less than infinity, it is used to get an approximation
for the scale of the argument to the function.

$head x_upper$$
This is a lower limit for the argument $icode x_in$$ to $icode fun$$.
If it is greater than - infinity, it is used to get an approximation
for the scale of the argument to the function.

$head x$$
Is the argument value at which the derivative is checked.

$head f$$
Is the value of the function at $icode x$$.

$head dfdx$$
Is the value of the derivative that we are checking.

$head rel_tol$$
Is an acceptable relative tolerance for the difference between
$icode dfdx$$ and its finite difference approximation.

$head result$$
The structure $icode%result%[%i%]%$$ contains
the following information about the $th i$$ component of the function:

$subhead rel_err$$
The smallest relative error that $code one_dim_derivative_check$$ found,
for the $th i$$ component of the function.
If it is less than or equal $icode rel_tol$$ the derivative check passed.
Searching for a step size with a smaller relative error stops as soon
as the derivative check passes for all components of the function.
If $icode%result%[%i%].rel_err%$$ is infinity,
then we were not able to evaluate the function.

$subhead step$$
The step size corresponding to $icode%result%[%i%].rel_err%$$.

$subhead dfdx$$
The finite difference approximation for the derivative,
corresponding to $icode%result%[%i%].rel_err%$$.

$end
*/
# include <cppad/mixed/typedef.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
// BEGIN_PROTOTYPE
struct one_dim_derivative_result {
	double rel_err;
	double step;
	double apx_dfdx;
};
template <class Object>
CppAD::vector<one_dim_derivative_result> one_dim_derivative_chk(
	Object&    obj          ,
	double     x_lower      ,
	double     x_upper      ,
	double     x            ,
	d_vector   f            ,
	d_vector   dfdx         ,
	double     rel_tol      )
// END_PROTOTYPE
{	// infinity
	double infinity = std::numeric_limits<double>::infinity();
	//
	// nan
	double nan      = std::numeric_limits<double>::quiet_NaN();
	//
	// sqrt_eps
	double sqrt_eps = std::sqrt( std::numeric_limits<double>::epsilon() );
	//
	// m
	size_t m = f.size();
	assert( m == dfdx.size() );
	//
	// x_max_abs
	double x_max_abs = std::fabs(x);
	if( std::fabs(x_upper) < infinity )
		x_max_abs = std::max(x_max_abs, std::fabs(x_upper));
	if( std::fabs(x_lower) < infinity )
		x_max_abs = std::max(x_max_abs, std::fabs(x_lower));
	//
	// log_max_rel_step
	// log of the maximum relative step to try
	double log_max_rel_step = std::log(0.1);
	if( x_max_abs == 0.0 )
		log_max_rel_step = std::log(1e+5);
	//
	// log_min_rel_step
	// log of the minimum relative step to try
	double log_min_rel_step = std::log(1e-10);
	//
	// n_try
	// number of relative steps to try
	size_t n_try = 5;
	if( x_max_abs == 0.0 )
		n_try = 7;
	//
	// log_diff
	// difference of log of relative step between trys
	double log_diff = (log_max_rel_step - log_min_rel_step)/double(n_try - 1);
	//
	// result
	// initialize
	CppAD::vector<one_dim_derivative_result> result(m);
	for(size_t i = 0; i < m; ++i)
	{	result[i].rel_err   = infinity;
		result[i].step      = nan;
		result[i].apx_dfdx  = nan;
	}
	//
	// loop over finite difference step sizes
	size_t i_try       = 0;        // index for this step size
	double rel_err_max = infinity; // maximum in result
	d_vector f_plus(m), f_minus(m);
	while( i_try < n_try && rel_err_max > rel_tol )
	{	// rel_step
		double log_try  = log_max_rel_step - log_diff * double(i_try);
		double rel_step = std::exp(log_try);
		//
		// step
		double step     = rel_step;
		if( x_max_abs != 0.0 )
			step = rel_step * x_max_abs;
		//
		// x_plus
		double x_plus = std::min(x + step, x_upper);
		//
		// x_minus
		double x_minus = std::max(x - step, x_lower);
		//
		// f_plus
		bool ok = obj.one_dim_function(x_plus, f_plus);
		//
		// f_minus
		ok = ok && obj.one_dim_function(x_minus, f_minus);
		//
		// rel_err_max
		rel_err_max = 0.0;
		for(size_t i = 0; i < m; ++i)
		{	// f_max_abs
			double f_max_abs = std::fabs(f[i]);
			f_max_abs        = std::max(f_max_abs, std::fabs(f_plus[i]));
			f_max_abs        = std::max(f_max_abs, std::fabs(f_minus[i]));
			//
			// apx_dfdx
			// finite difference approximation for derivative
			double apx_dfdx = (f_plus[i] - f_minus[i]) / (x_plus - x_minus);
			//
			// rel_err
			double rel_err = std::fabs(dfdx[i] - apx_dfdx);
			double den     = std::fabs(dfdx[i]) + std::fabs(apx_dfdx);
			den           += sqrt_eps * f_max_abs;
			if( den > 0.0 )
				rel_err = rel_err / den;
			//
			// best so far ?
			if( ok && 1.1 * rel_err < result[i].rel_err )
			{	result[i].rel_err  = rel_err;
				result[i].step     = step;
				result[i].apx_dfdx = apx_dfdx;
			}
			// rel_err_max
			rel_err_max = std::max(rel_err_max, result[i].rel_err);
		}
		//
		// next try
		++i_try;
	}
	return result;
}
} } // END_CPPAD_MIXED_NAMESPACE

# endif
