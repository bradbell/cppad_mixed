// $Id:$
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

$head Under Construction$$

$head Syntax$$
$icode%status% = CppAD::mixed::box_newton(
	%option%, %objective%, %x_low%, %x_up%, %x_in%, %x_out%
)%$$

$head Prototype$$
$srcfile%src/box_newton.cpp%0%// BEGIN PROTOTYPE%// END PROTOTYPE%1%$$

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
$srcfile%src/box_newton.cpp%0%// BEGIN OPTION%// END OPTION%1%$$

$subhead tolerance$$
This is the convergence tolerance for the optimization. The method has
converged when the absolute value of the
change in the current best approximate solution
is less than or equal $icode%option%.tolerance%$$ for all components
of $latex x$$.

$subhead print_level$$
This is the level of printing during this optimization process.
The default value $code 0$$ corresponds to no printing.

$subhead max_iter$$
This is the maximum number of iterations for the algorithm.
Each iterations of the algorithm corresponds to one call to
$codei%
	%w% = %objective%.solve(%x%, %v%)
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
	%v% = %solve%.solve(%x%, %w%)
%$$

$subhead x$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %x%
%$$
and is the same as in the previous call to
$icode%objective%.fun(%x%)%$$.

$subhead w$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %w%
%$$
and size $icode n$$.

$subhead v$$
The return value has prototype
$codei%
	CppAD::vector<double> %v%
%$$
It size is $code n$$ and it solves the equation
$latex \[
	f^{(2)} ( x ) \; v = w
\] $$
This equation can be solved because
we assume that $latex f^{(2)} (x)$$ is positive definite,

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
$srcfile%src/box_newton.cpp%0%// BEGIN STATUS%// END STATUS%1%$$


$end
------------------------------------------------------------------------------
*/
# include <Eigen/Core>
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
	// v^T * Hessian * v not positive where v corresponds to previous solve
	box_newton_neg_enum
};
// END STATUS

// BEGIN PROTOTYPE
template <class Objective>
box_newton_status box_newton(
	const box_newton_option&                        option     ,
	Objective&                                      objective  ,
	const Eigen::Matrix<double, Eigen::Dynamic, 1>& x_low      ,
	const Eigen::Matrix<double, Eigen::Dynamic, 1>& x_up       ,
	const Eigen::Matrix<double, Eigen::Dynamic, 1>& x_in       ,
	const Eigen::Matrix<double, Eigen::Dynamic, 1>& x_out      )
// END PROTOTYPE
{	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> double_vec;
	typedef Eigen::Matrix<bool, Eigen::Dynamic, 1>   bool_vec;

	size_t n = x_low.size();
	assert( n == x_up.size() );
	assert( n == x_in.size() );
	assert( n == x_out.size() );

	// initialize
	double_vec x_cur(n), g_cur(n), p_cur(n), d_cur(n), dx_cur(n), x_next(n);
	bool_vec active(n);
	double f_cur, eps;
	x_cur  = x_in;
	f_cur  = objective.fun(x_cur);
	g_cur  = objective.grad(x_cur);
	x_out  = x_cur;
	eps    = 100. * std::numeric_limits<double>::epsilon();
	//
	size_t iter = 0;
	while(iter < option.max_iter )
	{	iter++;
		//
		// set active set and projected gradient
		p_cur = g_cur;
		for(size_t i = 0; i < n; i++)
		{	bool lower = CppAD::NearEqual(x_cur[i], x_low[i], eps, eps);
			lower     &= g_cur[i] > 0.0;
			bool upper = CppAD::NearEqual(x_cur[i], x_up[i], eps, eps);
			lower     &= g_cur[i] < 0.0;
			active[i]  = lower || upper;
			if( active[i] )
				p_cur[i] = 0.0;
		}
		//
		// Netwon direction corresponding to projected gradient
		d_cur  = objective.solve(x_cur, p_cur);
		//
		// check for convergence
		double dx_norm = 0.0;
		for(size_t i = 0; i < n; i++)
		{	double xi = x_cur[i] - d_cur[i];
			xi        = std::max(xi, x_low[i]);
			xi        = std::min(xi, x_up[i]);
			dx_cur[i] = xi - x_cur[i];
			dx_norm   = std::max( dx_norm, std::fabs(dx_cur[i]) );
		}
		if( dx_norm < option.tolerance )
			return box_newton_ok_enum;
		//
		// directionanl derivative in dx_cur and p_cur directions
		double df_dx = g_cur.transpose() * dx_cur;
		double df_p  = g_cur.transpose() * p_cur;
		if( df_p <= 0.0 )
			return box_newton_ok_enum;
		//
		// if df_dx is not negative enough, use - p_cur direction
		if( df_dx > - df_p / 10. )
		{	dx_cur = - p_cur;
			df_dx  = - df_p;
		}
		//
		// line search
		double lam   = 1.0;
		size_t count = 0;
		double rate  = 0.0;
		while( count < option.max_line && rate > df_dx / 10. )
		{	count++;
			x_next = x_cur + lam * dx_cur;
			double f_next = objective.fun(x_next);
			double rate   = (f_next - f_cur) / lam;
			lam           = lam / 2.0;
		}
		if( rate > df_dx / 10. )
			return box_newton_max_line_enum;
		//
		x_out = x_next;
	}
	return box_newton_max_iter_enum;
}

} } // END_CPPAD_MIXED_NAMEPSPACE
