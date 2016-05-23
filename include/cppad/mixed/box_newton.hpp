// $Id:$
# ifndef CPPAD_MIXED_BOX_NEWTON_HPP
# define CPPAD_MIXED_BOX_NEWTON_HPP
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin box_newton$$
$spell
	enum
	CppAD
	const
	iter
$$

$section Newton's Optimization with Box Constraints$$

$head Syntax$$
$icode%status% = CppAD::mixed::box_newton(
	%option%, %objective%, %x_low%, %x_up%, %x_in%, %x_out%
)%$$

$head Prototype$$
$srcfile%include/cppad/mixed/box_newton.hpp
	%0%// BEGIN PROTOTYPE%// END PROTOTYPE%1%$$

$head Private$$
This function is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
Given a smooth function $latex f : \B{R}^n \rightarrow \B{R}$$,
together with its derivative and positive definite Hessian,
this routine solves the problem
$latex \[
	\R{minimize} \; f(x) \; \R{subject \; to} \; \ell \leq x \leq u
\] $$

$head option$$
$srcfile%include/cppad/mixed/box_newton.hpp
	%0%// BEGIN OPTION%// END OPTION%1%$$

$subhead tolerance$$
This is the convergence tolerance for the optimization. The method has
converged when the maximum absolute Newton step component
is less than or equal $icode%option%.tolerance%$$,
or when the derivative in the direction of the negative
projected gradient is not negative.

$subhead print_level$$
This is the level of printing during this optimization process.
The default value $code 0$$ corresponds to no printing.
$codei%
%print_level% >= 1
%$$
the iteration counter $icode iter$$,
the value of the objective $icode f$$,
maximum absolute component in the Newton step $icode |d|$$,
and line search step size $icode lam$$,
are  printed at the end of each iteration.
In addition, the return $icode status$$ is printed.
Note that a line search step size starts of $code 1$$
means that a Newton step was take for this iteration.
$codei%
%print_level% >= 2
%$$
the current argument vector $icode x$$ is printed for each iteration.
$codei%
%print_level% >= 3
%$$
the gradient $icode g$$ and
the negative of the gradient projected onto the constraint box
$icode p$$ are printed for each iteration.
$codei%
%print_level% >= 4
%$$
the Newton step $icode d$$ is printed for each iteration.

$subhead max_iter$$
This is the maximum number of iterations for the algorithm.
Each iterations of the algorithm corresponds to one call to
$codei%
	%d% = %objective%.solve(%x%, %p%)
%$$

$subhead max_line$$
This is the maximum number of line search steps to
perform during each iteration.
Each line search step corresponds to one call to
$codei%
	%f% = %objective%.fun(%x%)
%$$

$head fun$$
The object $icode objective$$ supports the following syntax
$codei%
	%f% = %objective%.fun(%x%)
%$$

$subhead x$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %x%
%$$
and size $icode n$$ and
is within the limits $icode x_low$$ and $icode x_up$$.
It specifies the argument value at which we are evaluating
the function.

$subhead f$$
The return value has prototype
$codei%
	double %f%
%$$
and is the value of $latex f(x)$$.

$head grad$$
The object $icode objective$$ supports the following syntax
$codei%
	%g% = %objective%.grad(%x%)
%$$

$subhead x$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %x%
%$$
and is the same as in the previous call to
$icode%objective%.fun(%x%)%$$.

$subhead g$$
The return value has prototype
$codei%
	CppAD::vector<double> %g%
%$$
Its size is $icode n$$
and it is the value of the gradient $latex f^{(1)} (x)^\R{T}$$.

$head solve$$
The object $icode objective$$ supports the following syntax
$codei%
	%d% = %solve%.solve(%x%, %p%)
%$$

$subhead x$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %x%
%$$
and is the same as in the previous call to
$icode%objective%.fun(%x%)%$$.

$subhead p$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %p%
%$$
and size $icode n$$.
This is the negative of the gradient projected onto the constraint box.

$subhead d$$
The return value has prototype
$codei%
	CppAD::vector<double> %d%
%$$
It size is $code n$$ and it solves the equation
$latex \[
	f^{(2)} ( x ) \; d = p
\] $$
This equation can be solved because
we assume that $latex f^{(2)} (x)$$ is positive definite,
The value $icode d$$ is the Newton step.

$head x_low$$
This vector has size $icode n$$ and specifies the lower limits
$latex \ell$$. If $icode%x_low%[%j%]%$$ is minus infinity,
there is no lower limit for the corresponding component of $latex x$$.

$head x_up$$
This vector has size $icode n$$ and specifies the upper limits
$latex u$$. If $icode%x_up%[%j%]%$$ is plus infinity,
there is no upper limit for the corresponding component of $latex x$$.

$head x_in$$
This vector has size $icode n$$ and specifies the initial value for
$latex u$$ during the optimization.
It must hold that, for $icode%j% = 0 , %...%, %n%-1%$$,
$codei%
	%x_low%[%j%] <= %x_in%[%j%] <= %x_up%[%j%]
%$$

$head x_out$$
This vector has size is $code n$$.
The input value of its elements does not matter.
Upon return it is the best approximate solution so far.
If $icode status$$ is $code box_newton_ok_enum$$,
for $latex j = 0, \ldots , n-1$$ the absolute value of the
change in $icode%x_out%[%j%]%$$ (between iterations)
is less than or equal $icode%option%.tolerance%$$.

$head status$$
The return value is one of the following enum values
$srcfile%include/cppad/mixed/box_newton.hpp
	%0%// BEGIN STATUS%// END STATUS%1%$$


$end
------------------------------------------------------------------------------
*/
# include <cppad/utility/vector.hpp>
# include <cppad/utility/near_equal.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN OPTION
struct box_newton_option {
	double tolerance;
	size_t print_level;
	size_t max_iter;
	size_t max_line;
	box_newton_option(void) : // set default values
	tolerance(1e-10)    ,
	print_level(0)      ,
	max_iter(50)        ,
	max_line(20)
	{}
};
// END OPTION

// BEGIN STATUS
enum box_newton_status {
	box_newton_ok_enum       , // x_out is ok
	box_newton_max_iter_enum , // maximum number of iterations reached
	box_newton_max_line_enum , // maximum number of line search steps reached
	// d^T * f^{(2)}(x) * d not positive where d is the Newton step
	box_newton_neg_enum
};
// END STATUS

// BEGIN PROTOTYPE
template <class Objective>
box_newton_status box_newton(
	const box_newton_option&      option     ,
	Objective&                    objective  ,
	const CppAD::vector<double>&  x_low      ,
	const CppAD::vector<double>&  x_up       ,
	const CppAD::vector<double>&  x_in       ,
	CppAD::vector<double>&        x_out      )
// END PROTOTYPE
{	size_t n = x_low.size();
	assert( n == x_up.size() );
	assert( n == x_in.size() );
	assert( n == x_out.size() );
# ifndef NDEBUG
	for(size_t i = 0; i < n; i++)
	{	assert( x_low[i] <= x_in[i] );
		assert( x_in[i]  <= x_up[i] );
	}
# endif
	// initialize
	CppAD::vector<double> x_cur(n), g_cur(n), p_cur(n), d_cur(n), dx_cur(n);
	CppAD::vector<double> x_next(n);
	CppAD::vector<bool> active(n);
	double f_cur, eps, inf;
	x_cur  = x_in;
	f_cur  = objective.fun(x_cur);
	eps    = 100. * std::numeric_limits<double>::epsilon();
	inf    = std::numeric_limits<double>::infinity();
	//
	size_t iter     = 0;
	while(iter < option.max_iter )
	{	iter++;
		if( option.print_level >= 2 )
			std::cout << "x = " << x_cur << std::endl;
		//
		// current gradient
		g_cur = objective.grad(x_cur);
		//
		// set active set and projected gradient
		for(size_t i = 0; i < n; i++)
		{	bool lower = false;
			if( x_low[i] > - inf )
			{	lower  = CppAD::NearEqual(x_cur[i], x_low[i], eps, eps);
				lower &= g_cur[i] > 0.0;
			}
			bool upper = false;
			if( x_up[i] < + inf )
			{	upper  = CppAD::NearEqual(x_cur[i], x_up[i], eps, eps);
				upper &= g_cur[i] < 0.0;
			}
			active[i]  = lower || upper;
			if( active[i] )
				p_cur[i] = 0.0;
			else
				p_cur[i] = - g_cur[i];
		}
		if( option.print_level >= 3 )
		{	std::cout << "g  = " << g_cur << std::endl;
			std::cout << "p  = " << p_cur << std::endl;
		}
		//
		// Netwon direction corresponding to projected gradient
		d_cur  = objective.solve(x_cur, p_cur);
		if( option.print_level >= 4 )
			std::cout << "d  = " << d_cur << std::endl;
		//
		// check for convergence
		double d_norm = 0.0;
		for(size_t i = 0; i < n; i++)
		{	double xi = x_cur[i] + d_cur[i];
			xi        = std::max(xi, x_low[i]);
			xi        = std::min(xi, x_up[i]);
			dx_cur[i] = xi - x_cur[i];
			d_norm    = std::max(d_norm, std::fabs(d_cur[i]) );
		}
		if( d_norm < option.tolerance )
		{	x_out = x_cur;
			std::cout << "box_newton_ok" << std::endl;
			return box_newton_ok_enum;
		}
		//
		// directionanl derivative in dx_cur and p_cur directions
		double df_dx = 0.0;
		double df_p  = 0.0;
		for(size_t i = 0; i < n; i++)
		{	df_dx += g_cur[i] * dx_cur[i];
			df_p  += g_cur[i] * p_cur[i];
		}
		if( df_p >= 0.0 )
		{	x_out = x_cur;
			std::cout << "box_newton_ok" << std::endl;
			return box_newton_ok_enum;
		}
		//
		// if df_dx is not negative enough, use - p_cur direction
		if( df_dx > df_p / 10. )
		{	for(size_t i = 0; i < n; i++)
				dx_cur[i] = p_cur[i];
			df_dx = df_p;
		}
		//
		// line search
		double lam   = 2.0;
		size_t count = 0;
		double rate  = 0.0;
		double f_next;
		while( count < option.max_line && rate > df_dx / 10. )
		{	count++;
			lam  = lam / 2.0;
			for(size_t i = 0; i < n; i++)
			{	x_next[i] = x_cur[i] + lam * dx_cur[i];
				x_next[i] = std::max(x_next[i], x_low[i]);
				x_next[i] = std::min(x_next[i], x_up[i]);
			}
			f_next  = objective.fun(x_next);
			rate    = (f_next - f_cur) / lam;
		}
		if( rate > df_dx / 10. )
		{	x_out = x_cur;
			std::cout << "box_newton_max_line" << std::endl;
			return box_newton_max_line_enum;
		}
		if( option.print_level > 1 )
			std::cout << std::endl;
		if( option.print_level >= 1 )
			std::cout << "iter = " << iter
			<< ", f = " << f_cur
			<< ", |d| = " << d_norm
			<< ", lam = " << lam << std::endl;
		//
		x_cur    = x_next;
		f_cur    = f_next;
	}
	x_out = x_cur;
	std::cout << "box_newton_max_iter" << std::endl;
	return box_newton_max_iter_enum;
}


} } // END_CPPAD_MIXED_NAMEPSPACE

# endif
