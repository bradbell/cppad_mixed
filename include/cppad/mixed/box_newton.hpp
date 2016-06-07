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
$srcfile%include/cppad/mixed/box_newton_options.hpp
	%0%// BEGIN OPTION%// END OPTION%0%$$

$subhead tolerance$$
This is the convergence tolerance for the optimization. The method has
converged when the norm of the projected gradient is less than
$icode%tolerance% > 0%$$.

$subhead direction_ratio$$
If the derivative of the function value in the Newton step direction,
divided by its derivative in the scaled negative projected gradient direction,
is less than $icode direction_ratio$$,
the scaled negative projected gradient is used for the line search direction.
If $latex p$$ is the projected gradient, and $latex d$$ is the Newton step
direction, $latex s = - p |d| / |p|$$
is the scaled negative projected gradient.
Note that $icode direction_ratio$$ should be
greater than zero and less than or equal one.
(If it is greater than one, the scaled negative projected gradient
direction will always be used).

$subhead line_ratio$$
The step size in the line search direction will be decreases until
the descent of the objective, divided by the line search step size,
is less than or equal $icode line_ratio$$
times the initial directional derivative in the line search direction.
Note that this ratio should be greater than zero and less than one half
(one half is the average descent rate for a convex quadratic function).
If the change in the function value is close to
numerical precision, the norm of the projected gradient,
divided by the line search step size,
is compared to $icode line_ratio$$ and used for line search termination.

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
the derivative in the current $icode dx/|dx|$$ direction $icode f_dx$$,
the derivative in the current $icode s/|s|$$ direction $icode f_s$$,
and line search step size $icode lam_dx$$ or $icode lam_s$$,
are printed at the end of each iteration.
If $icode lam_dx$$ ($icode lam_s$$) is printed,
$icode dx$$ ($icode s$$) is used for the search direction.
In addition, the return $icode status$$ is printed.
If the line search step size is $code 1$$,
and $code lam_dx$$ is printed,
a Newton step was taken for this iteration.
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
the negative projected gradient $icode p$$,
the Newton step $icode d$$,
and its projection to the feasible set $icode dx$$,
and the scaled negative projected gradient $icode s$$,
are printed for each iteration.
$codei%
%print_level% >= 3
%$$
the line search parameter $icode lam$$,
corresponding function values $icode f$$,
and the average derivative in the current line search direction $icode f_lam$$,
$latex \[
	f_\lambda = \frac{ f(x + \lambda q ) - f(x) } { \lambda | q |}
\] $$
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
This is the scaled negative gradient projected onto the constraint box.

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
condition has been satisfied.

$head status$$
The return value is one of the following enum values
(which are in the $code CppAD::mixed$$ namespace).
$srcfile%include/cppad/mixed/box_newton_status.hpp
	%0%// BEGIN STATUS%// END STATUS%0%$$

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
# include <cppad/mixed/box_newton_options.hpp>
# include <cppad/mixed/box_newton_status.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE


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
	CppAD::vector<double> s_cur(n), x_next(n), p_next(n);
	CppAD::vector<bool> active(n);
	double f_cur, eps, small, inf;
	x_cur  = x_in;
	f_cur  = objective.fun(x_cur);
	eps    = 1e2 * std::numeric_limits<double>::epsilon();
	small  = 1e2 * std::numeric_limits<double>::min();
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
		// set active set and negative projected gradient
		for(size_t i = 0; i < n; i++)
		{	bool lower = false;
			if( x_low[i] > - inf )
			{	lower  = CppAD::NearEqual(x_cur[i], x_low[i], eps, small);
				lower &= g_cur[i] > 0.0;
			}
			bool upper = false;
			if( x_up[i] < + inf )
			{	upper  = CppAD::NearEqual(x_cur[i], x_up[i], eps, small);
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
		//
		// directional derivatvies and norms
		double f_dx    = 0.0;
		double f_p     = 0.0;
		double d_norm  = 0.0;
		double dx_norm = 0.0;
		double p_norm  = 0.0;
		for(size_t i = 0; i < n; i++)
		{	f_dx     += g_cur[i] * dx_cur[i];
			f_p      += g_cur[i] * p_cur[i];
			d_norm   += d_cur[i] * d_cur[i];
			dx_norm  += dx_cur[i] * dx_cur[i];
			p_norm   += p_cur[i] * p_cur[i];
		}
		assert( f_p == - p_norm );
		d_norm  = std::sqrt( d_norm );
		dx_norm = std::sqrt( dx_norm );
		p_norm  = std::sqrt( p_norm );
		if( dx_norm > 0.0 )
			f_dx = f_dx / dx_norm;
		if( p_norm > 0.0 )
			f_p  = f_p / p_norm;
		//
		// check for convergence
		if( p_norm < options.tolerance )
		{	x_out = x_cur;
			if( options.print_level >= 1 )
				std::cout << "box_newton_ok: |p| = " << p_norm << std::endl;
			return box_newton_ok_enum;
		}
		//
		// scaled negative projected gradient
		for(size_t i = 0; i < n; i++)
			s_cur[i] = p_cur[i] * d_norm / p_norm;
		double f_s    = f_p;
		double s_norm = d_norm;
		//
		if( options.print_level >= 2 )
		{	std::cout << "x  = " << x_cur << std::endl;
			std::cout << "g  = " << g_cur << std::endl;
			std::cout << "p  = " << p_cur << std::endl;
			std::cout << "d  = " << d_cur << std::endl;
			std::cout << "dx = " << dx_cur << std::endl;
			std::cout << "s  = " << s_cur << std::endl;
		}
		//
		// if f_dx is not negative enough, use s_cur direction
		bool use_s    = f_dx > f_s * options.direction_ratio;
		double f_q    = f_dx;
		double q_norm = dx_norm;
		if( use_s )
		{	f_q    = f_s;
			q_norm = s_norm;
		}
		//
		// line search
		double lam   = 4.0;
		size_t count = 0;
		double f_lam = 0.0;
		double p_lam = 0.0;
		double f_next;
		while(
			count < options.max_line               &&
			f_lam > f_q * options.line_ratio       &&
			p_lam > - p_norm * options.line_ratio  )
		{	count++;
			lam  = lam / 4.0;
			for(size_t i = 0; i < n; i++)
			{	if( use_s )
					x_next[i] = x_cur[i] + lam * s_cur[i];
				else
					x_next[i] = x_cur[i] + lam * dx_cur[i];
				x_next[i] = std::max(x_next[i], x_low[i]);
				x_next[i] = std::min(x_next[i], x_up[i]);
			}
			f_next  = objective.fun(x_next);
			f_lam   = (f_next - f_cur) / (lam * q_norm);
			p_lam   = 0.0;
			if( CppAD::NearEqual(f_cur, f_next, eps, small) )
			{	// gradient
				p_next  = objective.grad(x_next);
				for(size_t i = 0; i < n; i++)
				{	// Quick way to convert gradient to projected gradient.
					// Not quiet correct and perhaps we should compute
					// projection for the point x_next.
					if( p_cur[i] == 0.0 )
						p_next[i] = 0.0;
					// norm squared of projected gradient
					p_lam += p_next[i] * p_next[i];
				}
				// norm of projected gradient at x_next
				p_lam = std::sqrt( p_lam );
				// rate of descent of projected gradient
				p_lam = (p_lam - p_norm) / (lam * q_norm);
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
				<< ", f_s = " << f_s;
			if( use_s )
				std::cout << ", lam_s = "  << lam << std::endl;
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
