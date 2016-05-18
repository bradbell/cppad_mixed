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
	CppAD
	const
	iter
$$

$section Newton's Optimization with Box Constraints$$

$head Under Construction$$

$head Syntax$$
$icode%x_out% = CppAD::mixed::box_newton(
	%x_low%, %x_up%, %x_in%, %options%, %objective%
)%$$

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

$head Prototype$$
$srcfile%src/box_newton.cpp%0%// BEGIN PROTOTYPE%// END PROTOTYPE%1%$$

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

$head options$$
This argument has prototype
$codei%
	const CppAD::mixed::box_newton_options& %options%
%$$

$subhead tolerance$$
This option has prototype
$codei%
	double %options%.tolerance
%$$
It is the convergence tolerance for the optimization. The method has
converged when the infinity norm of the projected gradient
$latex | p ( x ) |_\infty$$ is less than $icode%options%.tolerance%$$.
The default value for this option is $code 1e-10$$.

$subhead print_level$$
This option has prototype
$codei%
	size_t %options%.print_level
%$$
It is the level of printing during this optimization process.
The default value for this option is $code 0$$ which corresponds
to no printing.

$subhead max_iter$$
This option has prototype
$codei%
	size_t %options%.max_iter
%$$
It is the maximum number of iterations for the algorithm.
Each iterations of the algorithm corresponds to one call to
$codei%
	%w% = %objective%.solve(%x%, %v%)
%$$
The default value for this option is $code 50$$

$subhead max_line_search$$
This option has prototype
$codei%
	size_t %options%.max_line_search
%$$
It is the maximum number of line search steps to
perform during each iteration.
Each line search step corresponds to one call to
$codei%
	%f% = %objective%.fun(%x%)
%$$
The default value for this option is $code 20$$

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
and it is the value of the gradient of $latex f(x)$$.

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


$head x_out$$
The return value has prototype
$codei%
	CppAD::vector<double> %x_out%
%$$
It size is $code n$$ and it is a value of $icode x$$ such that
the infinity norm of the projected gradient
$latex | p ( x ) |_\infty$$ is less than $icode%options%.tolerance%$$.

$end
------------------------------------------------------------------------------
*/
namespace CppAD < namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

struct box_newton_options {
	double tolerance;
	size_t print_level;
	size_t max_iterions;
	size_t max_line_search;
	box_newton_options(void)
	: tolerance(1e-10)    ,
	: print_level(0)      ,
	: max_iter(50)        ,
	: max_line_search(20)
	{}
}

// BEGIN PROTOTYPE
template <class Objective>
CppAD::vector<double> box_newton(
	const CppAD::vector<double>&  x_low      ,
	const CppAD::vector<double>&  x_up       ,
	const CppAD::vector<double>&  x_in       ,
	box_newton_options            options    ,
	Objective&                    objective  )
// END PROTOTYPE




} } // END_CPPAD_MIXED_NAMEPSPACE
