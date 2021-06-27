/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin init_ldlt_ran_hes$$
$spell
	uu
	rc
	rcv
	ldlt_ran_hes
	init
	cppad
	CppAD
	eigen
	Cholesky
$$

$section Initialize Cholesky Factor of Hessian of Random Likelihood$$

$head Syntax$$
$icode%mixed_object%.init_ldlt_ran_hes()%$$

$head Private$$
This $code cppad_mixed$$ is a $cref private_base_class$$ member function.

$head Assumptions$$
The member variable
$cref/init_ran_hes_done_/init_ran_hes/init_ran_hes_done_/$$ is true.

$head init_ldlt_ran_hes_done_$$
The input value of this member variable must be false.
Upon return it is true.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head ldlt_ran_hes_$$
The member variable
$codei%
	CPPAD_MIXED_LDLT_CLASS ldlt_ran_hes_
%$$
must not been previously initialized
$codei%
	size_t ldlt_ran_hes_
%$$
Upon return, the function
$codei%
	ldlt_ran_hes_.init(ran_hes_uu_rcv_.pat())
%$$
has been called with $icode ran_hes_uu_rcv_.pat()$$
equal to the sparsity pattern for the
Hessian of the random likelihood with respect to the random effects;
see $cref ldlt_eigen_init$$.
This sparsity information is relative to just the random effects;
i.e., the row and column indices are relative to the vector $latex u )$$.
For each row and column indices in $icode ran_hes_uu_rcv_.pat()$$,
the corresponding row and column index in $code ran_hes_uu_rcv_$$ is
$code n_fixed_$$ greater.


$head a1_ldlt_ran_hes_$$
The member variable
$codei%
	CppAD::mixed::ldlt_eigen<a1_double> a1_ldlt_ran_hes_
%$$
is initialized the same as $code ldlt_ran_hes_$$.

$end
*/
# include <cppad/mixed/cppad_mixed.hpp>

void cppad_mixed::init_ldlt_ran_hes(void)
{	assert( ! init_ldlt_ran_hes_done_ );
	assert( init_ran_hes_done_ );
	//
	//
	// initialize ldlt_ran_hes_
	ldlt_ran_hes_.init(ran_hes_uu_rcv_.pat());
	//
	// initialize a1_ldlt_ran_hes_
	a1_ldlt_ran_hes_.init(ran_hes_uu_rcv_.pat());
	//
	init_ldlt_ran_hes_done_ = true;
	return;
}
