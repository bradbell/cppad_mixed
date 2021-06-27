/*
-----------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-----------------------------------------------------------------------------
$begin information_mat$$
$spell
	CppAD
	cppad
	rcv
$$

$section Compute the Observed Information For Fixed Effects$$

$head Deprecated 2020-03-22$$
Use $cref hes_fixed_obj$$ instead.

$head Syntax$$
$icode%information_rcv% = %mixed_object%.information_mat(
	%solution%, %random_opt%
)%$$

$head Purpose$$
Compute the observed information matrix.
We use $latex L ( \theta )$$ to denote the
$cref/fixed effects objective/theory/Objective/Fixed Effects Objective, L(theta)/$$.
The observed information is
$latex \[
	L^{(2)} ( \hat{\theta} )
\]$$
Absolute value terms in the
$cref/negative log-density vector/cppad_mixed/Negative Log-Density Vector/$$
for the $cref fix_likelihood$$ are not include in this Hessian
(because they do not have a derivative, let alone Hessian, at zero).

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head solution$$
is the $cref/solution/optimize_fixed/solution/$$
for a previous call to $cref optimize_fixed$$.
Only the $icode%solution%.fixed_opt%$$ field is used.

$head random_opt$$
is the optimal random effects corresponding to the solution; i.e.
$codei%
	%random_opt% = %mixed_object%.optimize_random(
		%random_options%,
		%solution%.fixed_opt,
		%random_lower%,
		%random_upper%,
		%random_in%
	)
%$$
$icode random_options$$,
$icode random_lower$$,
$icode random_upper$$, and
$icode random_in$$, are the same
as in the call to $code optimize_fixed$$ that corresponds to $icode solution$$.

$head information_rcv$$
The return value has prototype
$codei%
	CppAD::mixed::d_sparse_rcv %information_rcv%
%$$
see $cref/d_sparse_rcv/typedef/Sparse Types/d_sparse_rcv/$$.
This is a sparse matrix representation for the
lower triangle of the observed information matrix,
which is symmetric and hence determined by its lower triangle.
Absolute value terms in the
$cref/negative log-density vector/cppad_mixed/Negative Log-Density Vector/$$
for the $cref fix_likelihood$$ are not include in this Hessian
because they do not have a derivative (let alone Hessian) at zero.

$children%
	example/user/information_mat.cpp
%$$

$head Example$$
The file $cref information_mat.cpp$$ contains an example and
test of this routine. It returns true for success and false for failure.

$end
------------------------------------------------------------------------------
*/
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/triple2eigen.hpp>
# include <cppad/mixed/exception.hpp>

// ---------------------------------------------------------------------------
CppAD::mixed::d_sparse_rcv cppad_mixed::information_mat(
	const CppAD::mixed::fixed_solution&  solution             ,
	const d_vector&                      random_opt           )
{	d_sparse_rcv result;
	try
	{	result = try_hes_fixed_obj(solution.fixed_opt, random_opt);
	}
	catch(const std::exception& e)
	{	std::string error_message = "information_mat: std::exception: ";
		error_message += e.what();
		fatal_error(error_message);
		assert(false);
	}
	catch(const CppAD::mixed::exception& e)
	{	std::string error_message = e.message("information_mat");
		fatal_error(error_message);
		assert(false);
	}
	//
	return result;
}
