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
	namespace
	xam
	optimizer
	enum
	CppAD
	const
	iter
$$

$section Newton's Optimization with Box Constraints$$

$head Syntax$$
$icode%status% = CppAD::mixed::box_newton(
	%options%, %objective%, %x_low%, %x_up%, %x_in%, %x_out%
)%$$

$head Prototype$$
$srcfile%include/cppad/mixed/box_newton.hpp
	%0%// BEGIN PROTOTYPE%// END PROTOTYPE%1%$$

$head Public$$
This function is part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
Given a smooth function $latex f : \B{R}^n \rightarrow \B{R}$$,
together with its derivative and Hessian,
this routine solves the problem
$latex \[
	\R{minimize} \; f(x) \; \R{subject \; to} \; \ell \leq x \leq u
\] $$

$head options$$
The following declaration is made in the $code CppAD::mixed$$ namespace:
$srcfile%include/cppad/mixed/box_newton.hpp
	%0%// BEGIN OPTION%// END OPTION%1%$$

$subhead tolerance$$
This is the convergence tolerance for the optimization. The method has
converged when the norm of the projected gradient is less than
$icode%tolerance% > 0%$$.

$subhead direction_ratio$$
If the derivative of the function value in the Newton step direction,
divided by its derivative in the negative projected gradient direction,
is less than $icode direction_ratio$$,
the negative projected gradient is used for the line search direction.
Note that the directional derivatives are normalized before making this
comparison; i.e., divided by the norm of the corresponding direction.
Also note that this ratio should be greater than zero and less than one
(If it is greater than one, the negative projected gradient
direction will always be used).

$subhead line_ratio$$
The step size in the line search direction will be decreases until
the descent of the objective, divided by the line search step size,
is less than or equal $icode line_ratio$$
times the initial directional derivative in the line search direction.
Note this direction derivative is not normalized before making the comparison.
Also note that this ratio should be greater than zero and less than one half
(one half is the average descent rate for a convex quadratic function).
If the change in the function value is close to
numerical precision, the norm of the projected gradient,
divided by the line search step size,
is compared to $icode line_ratio$$ and used for line search
termination.

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

$subhead print_level$$
This is the level of printing during this optimization process.
The default value $code 0$$ corresponds to no printing.
$codei%
%print_level% >= 1
%$$
the iteration counter $icode iter$$,
the value of the objective $icode f$$,
norm of the Newton step projected to feasible set $icode |dx|$$,
the directive in the current $icode dx$$ direction $icode f_dx$$,
the directive in the current $icode p$$ direction $icode f_p$$,
and line search step size $icode lam_dx$$ or $icode lam_p$$,
are  printed at the end of each iteration.
If $icode lam_dx$$ ($icode lam_p$$) is printed,
$icode dx$$ ($icode p$$) is used for the search direction.
In addition, the return $icode status$$ is printed.
Note that a line search step size starts of $code 1$$
means that a Newton step was take for this iteration.
If the status is $code box_newton_ok$$,
the convergence test that passed is also printed.
$codei%
%print_level% >= 2
%$$
The lower limit $icode x_low$$, upper limit $icode x_up$$,
initial point $icode x_in$$,
and initial function value $icode f_in$$,
are printed before the first iteration.
The current argument vector $icode x$$,
the gradient $icode g$$,
the negative of the gradient projected onto the constraint box $icode p$$,
the Newton step $icode d$$,
and its projection to the feasible set $icode dx$$,
are printed for each iteration.
$codei%
%print_level% >= 3
%$$
the line search parameter $icode lam$$,
corresponding function values $icode f$$,
and the average derivative in the current line search direction $icode f_lam$$,
are printed for each iteration of the line search.
If the change in the objective is to small (near numerical precision),
the rate of descent of the norm of the
projected gradient $icode p_lam$$ is printed
instead of $icode f_lam$$.

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
	H \; d = p
\] $$
where $latex H$$ is an approximation for the Hessian
$latex f^{(2)} (x)$$.
For example, one possible choice for $latex H$$ is the value of the
Hessian at the initial point $icode x_in$$.
Another possible choice for $latex H$$ is actual Hessian $latex f^{(2)} (x)$$.

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
the $cref/tolerance/box_newton/options/tolerance/$$
condition has been satisfied or the derivative in the
direction of the negative projected gradient is non-negative.

$head status$$
The return value is one of the following enum values
(which are in the $code CppAD::mixed$$ namespace).
$srcfile%include/cppad/mixed/box_newton.hpp
	%0%// BEGIN STATUS%// END STATUS%1%$$

$children%example/user/box_newton_xam.cpp
%$$
$head Example$$
The file $cref box_newton_xam.cpp$$ contains an example and test
of this optimizer.


$end
------------------------------------------------------------------------------
*/
# include <cppad/utility/vector.hpp>
# include <cppad/utility/near_equal.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN OPTION
struct box_newton_options {
	double tolerance;
	double direction_ratio;
	double line_ratio;
	size_t max_iter;
	size_t max_line;
	size_t print_level;
	box_newton_options(void) : // set default values
	tolerance(1e-6)       ,
	direction_ratio(0.1)  ,
	line_ratio(0.05)      ,
	max_iter(50)          ,
	max_line(10)          ,
	print_level(0)
	{}
};
// END OPTION

// BEGIN STATUS
enum box_newton_status {
	box_newton_ok_enum       , // x_out is ok
	box_newton_max_iter_enum , // maximum number of iterations reached
	box_newton_max_line_enum   // maximum number of line search steps reached
};
// END STATUS

// BEGIN PROTOTYPE
template <class Objective>
box_newton_status box_newton(
	const box_newton_options&      options   ,
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
	CppAD::vector<double> x_next(n), p_next(n);
	CppAD::vector<bool> active(n);
	double f_cur, eps, inf;
	x_cur  = x_in;
	f_cur  = objective.fun(x_cur);
	eps    = 1e2 * std::numeric_limits<double>::epsilon();
	inf    = std::numeric_limits<double>::infinity();
	//
	if( options.print_level >= 2 )
	{	std::cout << "x_low = " << x_low << std::endl;
		std::cout << "x_up  = " << x_up << std::endl;
		std::cout << "x_in  = " << x_in << std::endl;
		std::cout << "f_in  = " << f_cur << std::endl;
	}
	//
	size_t iter     = 0;
	while(iter < options.max_iter )
	{	iter++;
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
		//
		// Netwon direction corresponding to projected gradient
		d_cur  = objective.solve(x_cur, p_cur);
		//
		// dx is projection of d onto the constraints
		for(size_t i = 0; i < n; i++)
		{	double xi = x_cur[i] + d_cur[i];
			xi        = std::max(xi, x_low[i]);
			xi        = std::min(xi, x_up[i]);
			dx_cur[i] = xi - x_cur[i];
		}
		if( options.print_level >= 2 )
		{	std::cout << "x  = " << x_cur << std::endl;
			std::cout << "g  = " << g_cur << std::endl;
			std::cout << "p  = " << p_cur << std::endl;
			std::cout << "d  = " << d_cur << std::endl;
			std::cout << "dx = " << dx_cur << std::endl;
		}
		//
		// directional derivatvies and norms
		double f_dx    = 0.0;
		double f_p     = 0.0;
		double dx_norm = 0.0;
		double p_norm  = 0.0;
		for(size_t i = 0; i < n; i++)
		{	f_dx     += g_cur[i] * dx_cur[i];
			f_p      += g_cur[i] * p_cur[i];
			dx_norm  += dx_cur[i] * dx_cur[i];
			p_norm   += p_cur[i] * p_cur[i];
		}
		assert( f_p == - p_norm );
		dx_norm = std::sqrt( dx_norm );
		p_norm  = std::sqrt( p_norm );
		if( p_norm < options.tolerance )
		{	x_out = x_cur;
			if( options.print_level >= 1 )
				std::cout << "box_newton_ok: |p| = " << p_norm << std::endl;
			return box_newton_ok_enum;
		}
		//
		// if f_dx is not negative enough, use p_cur direction
		bool use_p = f_dx * p_norm / dx_norm > f_p * options.direction_ratio;
		double f_q = f_dx;
		if( use_p )
			f_q = f_p;
		//
		// line search
		double lam   = 2.0;
		size_t count = 0;
		double f_lam = 0.0;
		double p_lam = 0.0;
		double f_next;
		while(
			count < options.max_line               &&
			f_lam > f_q * options.line_ratio       &&
			p_lam > - p_norm * options.line_ratio  )
		{	count++;
			lam  = lam / 2.0;
			for(size_t i = 0; i < n; i++)
			{	if( use_p )
					x_next[i] = x_cur[i] + lam * p_cur[i];
				else
					x_next[i] = x_cur[i] + lam * dx_cur[i];
				x_next[i] = std::max(x_next[i], x_low[i]);
				x_next[i] = std::min(x_next[i], x_up[i]);
			}
			f_next  = objective.fun(x_next);
			f_lam   = (f_next - f_cur) / lam;
			p_lam   = 0.0;
			if( CppAD::NearEqual(f_cur, f_next, eps, eps) )
			{	// gradient
				p_next  = objective.grad(x_next);
				for(size_t i = 0; i < n; i++)
				{	// Quick way to convert gradient to projected gradient.
					// Note quiet correct and perhaps we should compute
					// projection for the point x_next.
					if( p_cur[i] == 0.0 )
						p_next[i] = 0.0;
					// norm squared of projected gradient
					p_lam += p_next[i] * p_next[i];
				}
				// norm of projected gradient at x_next
				p_lam = std::sqrt( p_lam );
				// rate of descent of projected gradient
				p_lam = (p_lam - p_norm) / lam;
				if( options.print_level >= 3 )
					std::cout << "lam = " << lam
					<< ", |p| = " << p_norm
					<< ", p_lam = " << p_lam
					<< std::endl;
			}
			else if( options.print_level >= 3 )
				std::cout << "lam = " << lam
				<< ", f = " << f_next
				<< ", f_lam = " << f_lam
				<< std::endl;
		}
		if( f_lam > f_q * options.line_ratio      &&
		    p_lam > - p_norm * options.line_ratio )
		{	x_out = x_cur;
			if( options.print_level >= 1 )
				std::cout << "box_newton_max_line" << std::endl;
			return box_newton_max_line_enum;
		}
		x_cur    = x_next;
		f_cur    = f_next;
		//
		if( options.print_level >= 1 )
		{	std::cout
				<< "iter = " << iter
				<< ", f = "    << f_cur
				<< ", |dx| = " << dx_norm
				<< ", f_dx = " << f_dx
				<< ", f_p = " << f_p;
			if( use_p )
				std::cout << ", lam_p = "  << lam << std::endl;
			else
				std::cout << ", lam_dx = "  << lam << std::endl;
		}
	}
	x_out = x_cur;
	if( options.print_level >= 1 )
			std::cout << "box_newton_max_iter" << std::endl;
	return box_newton_max_iter_enum;
}


} } // END_CPPAD_MIXED_NAMEPSPACE

# endif
