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

$head Purpose$$
Given a smooth function $latex f : \B{R}^n \rightarrow \B{R}$$,
together with its derivative and positive definite Hessian,
this routine solves the problem
$latex \[
	\R{minimize} \; f(x) \; \R{subject \; to} \; \ell \leq x \leq u
\] $$

$head Project Gradient$$
We define the projected gradient $latex p(x)$$ by
$latex \[
p_j (x) = \left\{ \begin{array}{ll}
	0 & \R{if} \; f_j^{(1)} ( x ) \geq 0 \; \R{and} \; x_j = \ell_j \\
	0 & \R{if} \; f_j^{(1)} ( x ) \leq 0 \; \R{and} \; x_j = u_j    \\
	f_j^{(1)} (x) & \R{otherwise}
	\end{array} \right.
\] $$
where $latex f_j^{(1)} (x)$$ is the $th j$$ component of the
derivative of $latex f(x)$$.

$head option$$
$srcfile%src/box_newton.cpp%0%// BEGIN OPTION%// END OPTION%1%$$

$subhead tolerance$$
This is the convergence tolerance for the optimization. The method has
converged when $latex | p_j ( x ) |$$ is less than
or equal $icode%option%.tolerance%$$ for $latex j = 0 , \ldots , n-1$$.

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
Upon return it is a value of $icode x$$ such that
$latex | p_j ( x ) |_\infty$$ is less than or equal
$icode%option%.tolerance%$$ for $latex j = 0, \ldots , n-1$$.

$head status$$
The return value is one of the following enum values
$srcfile%src/box_newton.cpp%0%// BEGIN STATUS%// END STATUS%1%$$


$end
------------------------------------------------------------------------------
*/
namespace CppAD < namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN OPTION
struct box_newton_option {
	double tolerance;
	size_t print_level;
	size_t max_iter;
	size_t max_line;
	box_newton_option(void) // set default values
	: tolerance(1e-10)    ,
	: print_level(0)      ,
	: max_iter(50)        ,
	: max_line(20)
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
	box_newton_option             option     ,
	Objective&                    objective  ,
	const CppAD::vector<double>&  x_low      ,
	const CppAD::vector<double>&  x_up       ,
	const CppAD::vector<double>&  x_in       ,
	CppAD::vector<double>&        x_out      )
// END PROTOTYPE




} } // END_CPPAD_MIXED_NAMEPSPACE
