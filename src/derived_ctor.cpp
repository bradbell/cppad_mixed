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
$begin derived_ctor$$
$spell
	CppAD
	bool
	cppad
	dismod
	const
	vec
	rcv
$$

$section User Defined Class Derived From cppad_mixed$$

$head Syntax$$
$icode%mixed_derived% %mixed_object%(
	%n_fixed%, %n_random%, %quasi_fixed%, %bool_sparsity%, %A_rcv%, %...%
)%$$

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

$head n_fixed$$
This argument has prototype
$codei%
	size_t %n_fixed%
%$$
and is the number of
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$ in the model.

$head n_random$$
This argument has prototype
$codei%
	size_t %n_random%
%$$
and is the number of
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$ in the model.
In the case where there are
$cref/no random effects/cppad_mixed/Problem/No Random Effects/$$,
$icode%n_random% = 0%$$.

$head quasi_fixed$$
This argument has prototype
$codei%
	bool %quasi_fixed%
%$$
Currently, true is the recommended value for this parameter
(expected to be faster and use less memory).
If $icode quasi_fixed$$ is true,
the a quasi-Newton approximation for the Hessian of the total objective
$cref/L(theta)/theory/Objective/Total Objective, L(theta)/$$
is used during the optimization of the fixed effects.
Otherwise, the Hessian of the total objective is computed using the
approximate Laplace objective
$cref/H(beta, theta, u)
	/theory
	/Approximate Laplace Objective, H(beta, theta, u)
/$$.
If $icode quasi_fixed$$ is true,
some initialization is skipped during $cref initialize$$ is skipped.
This initialization is needed, and hence done during
the computation of the $cref/information matrix/information_mat/$$.
If $icode quasi_fixed$$ is true,
the amount of memory used by the
$cref/mixed_derived/derived_ctor/mixed_derived/$$ object
after the information matrix is computed
is similar after then initialization when $icode quasi_fixed$$ is false.

$head bool_sparsity$$
This argument has prototype
$codei%
	bool %bool_sparsity%
%$$
If it is true, where possible
boolean sparsity patterns are used for this computation,
otherwise set sparsity patterns are used.
This should only affect to amount of time and memory used for the
computations.
If this argument is not present, the type of sparsity patterns
is not specified.

$head A_rcv$$
This argument has prototype
$codei%
	const CppAD::mixed::sparse_rcv& %A_rcv%
%$$
It is a
$cref/sparse matrix/sparse_mat_info/Notation/Sparse Matrix/$$
representation of the
$cref/random constraint matrix
	/cppad_mixed
	/Notation
	/Random Constraint Matrix, A
/$$
$latex A$$.
If $icode%random_vec%.size()%$$ is zero, this must be the empty matrix.
The member variable $icode A_rcv_$$ is set equal to $icode A_rcv$$
before any other routines are called by this routine.

$head ...$$
Other arguments to the derived class constructor
(that are not used by the base class constructor).
The other arguments need not appear at the end of the derived
class constructor (as in the syntax above).

$head CppAD ErrorHandler$$
If a CppAD error occurs, its
$href%http://www.coin-or.org/CppAD/Doc/errorhandler.htm%ErrorHandler%$$
is used to map it to either a
$cref/fatal_error/public/User Defined Functions/fatal_error/$$
or
$cref/warning/public/User Defined Functions/warning/$$.

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
	const CppAD::mixed::sparse_rcv&       A_rcv         )
:
n_fixed_(n_fixed)               ,
n_random_(n_random)             ,
quasi_fixed_(quasi_fixed)       ,
bool_sparsity_(bool_sparsity)   ,
A_rcv_(A_rcv)                 ,
init_ran_like_done_(false)      ,
init_ran_hes_done_(false)       ,
init_ldlt_ran_hes_done_(false)  ,
init_hes_cross_done_(false)     ,
init_newton_checkpoint_done_(false)   ,
init_laplace_obj_done_(false)    ,
init_laplace_obj_hes_done_(false),
init_fix_like_done_(false)      ,
init_fix_con_done_(false)       ,
initialize_done_(false)         ,
cppad_error_handler_(handler)   ,
ldlt_ran_hes_(n_random)
{ }

// base class destructor
cppad_mixed::~cppad_mixed(void)
{ }
