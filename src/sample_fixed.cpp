// $Id:$
/*
-----------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-----------------------------------------------------------------------------
$begin sample_fixed$$
$spell
	Cholesky
	covariance
	pairwise
	CppAD
$$

$section Simulating the Posterior Distribution for the Fixed Effects$$

$head Syntax$$
$icode%correlation% = %mixed_object%.sample_fixed(
	%sample%,
	%non_zero%,
	%solution%,
	%fixed_lower%,
	%fixed_upper%,
	%fixed_constraint_lower%,
	%fixed_constraint_upper%
)%$$

$head Under Construction$$

$head Purpose$$
Sample the asymptotic posterior distribution for the
optimal fixed effects (given the model and the data).

$head sample$$
This argument has prototype
$codei%
	CppAD::vector<double>& %sample%
%$$
and its size is
$cref/n_fixed/derived_ctor/n_fixed/$$ times $icode n_sample$$.
The input value of its elements does not matter.
Upon return,
for $codei%i% = 0 , %...%, %n_sample%-1%$$,
$codei%j% = 0 , %...%, %n_fixed%-1%$$,
$codei%
	%sample%[ %i% * %n_fixed% + %j% ]
%$$
is the $th j$$ component of the $th i$$ sample of the
optimal fixed effects $latex \hat{\theta}$$.

$head non_zero$$
This argument has prototype
$codei%
	double  %non_zero%
%$$
and is a value between zero and one.
It specifies the fraction of off diagonal elements in
the posterior covariance that are included in the simulation.

$head n_sample$$
This argument has prototype
$codei%
	size_t %n_sample%
%$$
and is the number of samples of the fixed effects to generate.

$head solution$$
is the $cref/solution/optimize_fixed/solution/$$
for a previous call to $cref optimize_fixed$$.

$head fixed_lower$$
is the $cref/fixed_lower/optimize_fixed/fixed_lower/$$
for a previous call to $code optimize_fixed$$.

$head fixed_upper$$
is the $cref/fixed_upper/optimize_fixed/fixed_upper/$$
for a previous call to $code optimize_fixed$$.

$head fix_constraint_lower$$
is the $cref/fix_constraint_lower/optimize_fixed/fix_constraint_lower/$$
for a previous call to $code optimize_fixed$$.

$head fix_constraint_upper$$
is the $cref/fix_constraint_upper/optimize_fixed/fix_constraint_upper/$$
for a previous call to $code optimize_fixed$$.

$head correlation$$
The return value has prototype
$codei%
	double %correlation%
%$$
This is the largest absolute correlation
that is treated as zero when simulating the posterior samples.
If $icode%non_zero% = 1.0%$$, no values are treated as zero
and $icode%correlation% = 0.0%$$.


$head Unconstrained Covariance$$
We use $latex L ( \theta )$$ to denote the
$cref/total objective/theory/Objective/Total Objective, L(theta)/$$.
If $latex \hat{\theta}$$ is the optimal value for the fixed effects,
the corresponding estimate for the covariance of $latex \hat{\theta}$$ is
given by
$latex \[
	\B{C} ( \hat{\theta}, \hat{\theta} )
	=
	L^{(2)} ( \hat{\theta} )^{-1}
\]$$.

$head Approximate Constraints$$
Let $latex n$$ be the number of fixed constraints,
$latex n$$ the number of active constraints,
and the equations $latex e( \theta ) = b$$ the active constraints
where $latex e : \B{R}^n \rightarrow \B{R}^m$$, $latex b \in \B{R}^m$$,
and the inequality constraints have been converted to equalities at the
active bounds.
Define the random variable $latex \tilde{e}$$ as the affine
approximation for $latex e( \theta )$$.
$latex \[
	\tilde{e} = b + e^{(1)} ( \hat{\theta} ) ( \theta - \hat{\theta} )
\] $$

$head Constrained Covariance$$
We approximate the distribution for $latex \hat{\theta}$$ as normal,
hence $latex \tilde{e}$$ is also normal.
Furthermore, we approximate the distribution for $latex \hat{\theta}$$,
by is conditional distribution given that $latex \tilde{e}$$
is equal to $latex b$$.
In other words,
$latex \[
	\B{C} ( \hat{\theta} | b \W{,} \hat{\theta} | b )
	=
	\B{C} ( \hat{\theta} \W{,} \hat{\theta} )
	-
	\B{C} ( \hat{\theta} \W{,} \tilde{e} )
	\B{C} ( \tilde{e}  \W{,} \tilde{e} )^{-1}
	\B{C} ( \tilde{e}  \W{,} \hat{\theta} )
\] $$
Using the notation
$latex C = L^{(2)} ( \hat{\theta} )^{-1}$$
and
$latex E = e^{(1)} ( \hat{\theta} )$$,
we define
$latex \[
	D
	=
	\B{C} ( \hat{\theta} | b \W{,} \hat{\theta} | b )
	=
	C - C E^\R{T} \left( E C E^\R{T} \right)^{-1}  E C
\] $$
We use $latex D( \alpha )$$ for the matrix where the
off diagonal elements of $latex D$$ corresponding to an absolute correlation
less than $latex \alpha$$ are replaced by zero.
To be specific,
$latex \[
D_{i,j} ( \alpha ) = \left\{ \begin{array}{ll}
0 & \R{if} \; i \neq j \; \R{and} \; \alpha \geq | D_{i,j} | / \sqrt{D_ii D_jj}
\\
D_{i,j}       & \R{otherwise}
\end{array} \right.
\] $$
We define $latex \bar{D} = D ( \alpha )$$ where $latex \alpha$$
is smallest value, greater than or equal zero,
such that the fraction of non-zero off diagonal values in
$latex \bar{D}$$ is greater than or equal $icode non_zero$$.
The return value $cref/correlation/sample_fixed/correlation/$$ is
the maximum value of $latex | D_{i,j} | / \sqrt{ D_ii D_jj }$$
that is converted to zero in the definition of $latex \bar{D}$$.

$head Simulation$$
We use $latex \sqrt{\bar{D}}$$ to denote a Cholesky factor of
$latex \bar{D}$$; i.e.,
$latex \bar{D} = \sqrt{\bar{D}} \sqrt{\bar{D}}^\R{T}$$.
Further suppose that for $icode%i% = 0 , %...%, %n_sample%-1%$$,
$latex w^i \in \B{R}^n$$ is simulated as independent normal random vectors
with identity covariance matrix.
The $th i$$ the samples for the fixed effects estimate $latex \hat{\theta}$$
is $latex \sqrt{\bar{D}} w^i$$.

$end
------------------------------------------------------------------------------
*/
# include <cppad/mixed/cppad_mixed.hpp>

double cppad_mixed::sample_fixed(
	d_vector&                            sample               ,
	double                               non_zero             ,
	const CppAD::mixed::fixed_solution&  solution             ,
	const d_vector&                      fixed_lower          ,
	const d_vector&                      fixed_upper          ,
	const d_vector&                      fix_constraint_lower ,
	const d_vector&                      fix_constraint_upper )
{	// under construction
	return 0.0;
}
