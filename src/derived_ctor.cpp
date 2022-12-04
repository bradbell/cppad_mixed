// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin derived_ctor}
{xrst_spell
   boolean
   nc
   nr
}

User Defined Class Derived From cppad_mixed
###########################################

Syntax
******

| *mixed_derived* *mixed_object* (
| |tab| *n_fixed* , *n_random* ,
| |tab| *quasi_fixed* , *bool_sparsity* , *A_rcv* , *trace_init* ,
| |tab| ...
| )

Prototype
*********
{xrst_literal
   include/cppad/mixed/base_class.hpp
   // BEGIN_CPPAD_MIXED_CTOR
   // END_CPPAD_MIXED_CTOR
}

See Also
********
:ref:`initialize-name`

mixed_derived
*************
This is the name of the class derived in the following fashion:

   ``class`` *mixed_derived* : ``public cppad_mixed`` {

mixed_object
************
This is the derived class object that is constructed by the syntax above.

cppad_mixed
***********
The derived class constructor must call its base class constructor as follows:

| |tab| ``cppad_mixed`` (
| |tab| |tab| *n_fixed* , *n_random* , *quasi_fixed* , *bool_sparsity* , *A_rcv*
| |tab| )

The arguments *quasi_fixed* , *bool_sparsity* , *A_rcv*
are optional; see default values in prototype above.

n_fixed
*******
This is the number of
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>` in the model.

n_random
********
This is the number of
:ref:`random effects<problem@Notation@Random Effects, u>` in the model.
In the case where there are
:ref:`problem@Maximum Likelihood@No Random Effects` ,
*n_random*  = 0 .

quasi_fixed
***********

true
====
If *quasi_fixed* is true,
a quasi-Newton approximation for the Hessian of the fixed effects objective
:ref:`L(theta)<theory@Objective@Fixed Effects Objective, L(theta)>`
is used during the optimization of the fixed effects.
This is more robust when :ref:`optimize_fixed@fixed_in`
is far away from a reasonable value and might lead to the
Hessian w.r.t. the random effects not being positive definite.
If *quasi_fixed* is true,
some initialization is skipped during :ref:`initialize-name` .
This initialization is needed, and hence computed if
the and when the :ref:`information matrix<information_mat-name>` is computed.

false
=====
If *quasi_fixed* is false,
the Hessian of the fixed effects objective is computed using the
approximate Laplace objective
:ref:`H(beta, theta, u)<theory@Approximate Laplace Objective, H(beta, theta, u)>` .
The extra routines for initializing the second order accurate
approximation for the Laplace objective ``init_laplace_obj_fun`` ,
and ``init_laplace_obj_hes`` are used to initialize the
Hessian of the fixed effects objective.

bool_sparsity
*************
If *bool_sparsity* is true, where possible
boolean sparsity patterns are used for this computation,
otherwise set sparsity patterns are used.
This should only affect to amount of time and memory used for the
computations.

A_rcv
*****
This is a
:ref:`sparse_mat_info@Notation@Sparse Matrix`
representation of the
:ref:`random constraint matrix<problem@Notation@Random Constraint Matrix, A>`
:math:`A`.
If *random_vec* . ``size`` () is zero,
there are no constraint equations and *A_rcv* . ``nr`` () == 0
(this is the case for the default value of this argument).
Otherwise, *A_rcv* . ``nc`` () must be equal to *n_random*
and *A_rcv* . ``nr`` () is the number of constraints.

trace_init
**********
If true, trace the initialization of cppad_mixed data structures on
standard output. This can be useful for large problems where the initialization
takes a significant amount of time.
For an example see
:ref:`hes_fixed_obj.cpp@trace_init` .

...
***
Other arguments to the derived class constructor
(that are not used by the base class constructor).
The other arguments need not appear at the end of the derived
class constructor (as in the syntax above).

CppAD ErrorHandler
******************
If a CppAD error occurs, its
`ErrorHandler <http://www.coin-or.org/CppAD/Doc/errorhandler.htm>`_
is used to map it to either a
:ref:`base_class@User Defined Functions@fatal_error`
or
:ref:`base_class@User Defined Functions@warning` .
{xrst_toc_hidden
   example/user/derived_ctor.cpp
}
Example
*******
The file :ref:`derived_ctor.cpp-name` contains an example and test
that uses this derived class.
It returns true for success and false for failure.

{xrst_end derived_ctor}
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
     {  // use the most recent cppad_mixed fatal_error routine
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
   const CppAD::mixed::d_sparse_rcv&     A_rcv         ,
   bool                                  trace_init    )
:
n_fixed_(n_fixed)                   ,
n_random_(n_random)                 ,
quasi_fixed_(quasi_fixed)           ,
bool_sparsity_(bool_sparsity)       ,
A_rcv_(A_rcv)                       ,
trace_init_(trace_init)             ,
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
{  if( A_rcv.nr() > 0 && A_rcv.nc() != n_random )
   {  std::string message = "cppad_mixed ctor: A_rcv.nr() > 0 and A_rcv.nc() != n_random";
      fatal_error(message);
   }
}

// base class destructor
cppad_mixed::~cppad_mixed(void)
{ }
