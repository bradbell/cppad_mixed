/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin derived_ctor$$
$spell
	CppAD
	bool
	cppad
	dismod
	const
	vec
	rcv
	init_obj_hes
	nr
	nc
$$

$section User Defined Class Derived From cppad_mixed$$

$head Syntax$$
$icode%mixed_derived% %mixed_object%(
	%n_fixed%, %n_random%, %quasi_fixed%, %bool_sparsity%, %A_rcv%, %...%
)%$$

$head Prototype$$
$srcfile%include/cppad/mixed/base_class.hpp%
	0%// BEGIN_CPPAD_MIXED_CTOR%// END_CPPAD_MIXED_CTOR%0
%$$

$head See Also$$
$cref initialize$$

$head mixed_derived$$
This is the name of the class derived in the following fashion:
$codei%
	class %mixed_derived% : public cppad_mixed {
%$$

$head mixed_object$$
This is the derived class object that is constructed by the syntax above.

$head cppad_mixed$$
The derived class constructor must call its base class constructor as follows:
$codei%
	cppad_mixed(
		%n_fixed%, %n_random%, %quasi_fixed%, %bool_sparsity%, %A_rcv%
	)
%$$
The arguments $icode quasi_fixed$$, $icode bool_sparsity$$, $icode A_rcv$$
are optional; see default values in prototype above.

$head n_fixed$$
This is the number of
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$ in the model.

$head n_random$$
This is the number of
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$ in the model.
In the case where there are
$cref/no random effects/cppad_mixed/Problem/No Random Effects/$$,
$icode%n_random% = 0%$$.

$head quasi_fixed$$

$subhead true$$
If $icode quasi_fixed$$ is true,
a quasi-Newton approximation for the Hessian of the fixed effects objective
$cref/L(theta)/theory/Objective/Fixed Effects Objective, L(theta)/$$
is used during the optimization of the fixed effects.
This is more robust when $cref/fixed_in/optimize_fixed/fixed_in/$$
is far away from a reasonable value and might lead to the
Hessian w.r.t. the random effects not being positive definite.
If $icode quasi_fixed$$ is true,
some initialization is skipped during $cref initialize$$.
This initialization is needed, and hence computed if
the and when the $cref/information matrix/information_mat/$$ is computed.

$subhead false$$
If $icode quasi_fixed$$ is false,
the Hessian of the fixed effects objective is computed using the
approximate Laplace objective
$cref/H(beta, theta, u)
	/theory
	/Approximate Laplace Objective, H(beta, theta, u)
/$$.
The extra routines for initializing the second order accurate
approximation for the Laplace objective $code init_laplace_obj_fun$$,
and $code init_laplace_obj_hes$$ are used to initialize the
Hessian of the fixed effects objective.

$head bool_sparsity$$
If $icode bool_sparsity$$ is true, where possible
boolean sparsity patterns are used for this computation,
otherwise set sparsity patterns are used.
This should only affect to amount of time and memory used for the
computations.

$head A_rcv$$
This is a
$cref/sparse matrix/sparse_mat_info/Notation/Sparse Matrix/$$
representation of the
$cref/random constraint matrix
	/cppad_mixed
	/Notation
	/Random Constraint Matrix, A
/$$
$latex A$$.
If $icode%random_vec%.size()%$$ is zero,
there are no constraint equations and $icode%A_rcv%.nr() == 0%$$
(this is the case for the default value of this argument).
Otherwise, $icode%A_rcv%.nc()%$$ must be equal to $icode n_random$$
and $icode%A_rcv%.nr()%$$ is the number of constraints.

$head ...$$
Other arguments to the derived class constructor
(that are not used by the base class constructor).
The other arguments need not appear at the end of the derived
class constructor (as in the syntax above).

$head CppAD ErrorHandler$$
If a CppAD error occurs, its
$href%http://www.coin-or.org/CppAD/Doc/errorhandler.htm%ErrorHandler%$$
is used to map it to either a
$cref/fatal_error/base_class/User Defined Functions/fatal_error/$$
or
$cref/warning/base_class/User Defined Functions/warning/$$.

$children%
	example/user/derived_ctor.cpp
%$$
$head Example$$
The file $cref derived_ctor.cpp$$ contains an example and test
that uses this derived class.
It returns true for success and false for failure.

$end
*/
# include<cppad/mixed/cppad_mixed.hpp>
# include<cppad/mixed/exception.hpp>

namespace { // BEGIN_EMPTY_NAMESPACE
	void handler(
		bool known       ,
		int  line        ,
		const char *file ,
		const char *exp  ,
		const char *msg  )
     {	// use the most recent cppad_mixed fatal_error routine
		std::string thrower, brief;
		//
		thrower += "\nCppAD: file = ";
		thrower += file;
		thrower += "\nline = ";
		thrower += CppAD::to_string(line);
		//
		brief   += "msg = ";
		brief   += msg;
		//
		// CppAD ErrorHandler: uncomment next line to debug CppAD asserts
		// assert(false);
		CppAD::mixed::exception e(thrower, brief);
		throw(e);
     }

} // END_EMPTY_NAMESPACE

// base class constructor
cppad_mixed::cppad_mixed(
	size_t                                n_fixed       ,
	size_t                                n_random      ,
	bool                                  quasi_fixed   ,
	bool                                  bool_sparsity ,
	const CppAD::mixed::d_sparse_rcv&     A_rcv         )
:
n_fixed_(n_fixed)                   ,
n_random_(n_random)                 ,
quasi_fixed_(quasi_fixed)           ,
bool_sparsity_(bool_sparsity)       ,
A_rcv_(A_rcv)                       ,
init_ran_like_done_(false)          ,
init_ran_jac_done_(false)           ,
init_ran_hes_done_(false)           ,
init_ldlt_ran_hes_done_(false)      ,
init_hes_cross_done_(false)         ,
init_laplace_obj_done_(false)       ,
init_laplace_obj_fun_done_(false)   ,
init_laplace_obj_hes_done_(false)   ,
init_fix_like_done_(false)          ,
init_fix_con_done_(false)           ,
initialize_done_(false)             ,
cppad_error_handler_(handler)       ,
ldlt_ran_hes_(n_random)             ,
a1_ldlt_ran_hes_(n_random)
{	if( A_rcv.nr() > 0 && A_rcv.nc() != n_random )
	{	std::string message = "cppad_mixed ctor: A_rcv.nr() > 0 and A_rcv.nc() != n_random";
		fatal_error(message);
	}
}

// base class destructor
cppad_mixed::~cppad_mixed(void)
{ }
