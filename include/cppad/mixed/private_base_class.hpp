// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_PRIVATE_BASE_CLASS_HPP
# define CPPAD_MIXED_PRIVATE_BASE_CLASS_HPP

// cppad_mixed.hpp: private member functions and variables
private:
/*
------------------------------------------------------------------------------
{xrst_begin private_base_class dev}
{xrst_spell
   boolean
   logdet
   uu
}

cppad_mixed: Private Declarations
#################################
These ``cppad_mixed`` class declarations are ``private`` .
They are **not** part of the user API,
they may change with time, and they
can **not** be used by a derived class object
:ref:`derived_ctor@mixed_object` .

Contents
********
{xrst_toc_table
   include/cppad/mixed/pack.hpp
   include/cppad/mixed/unpack.hpp
   src/init_ran_like.cpp
   src/init_ran_jac.cpp
   src/init_ran_hes.cpp
   src/init_ldlt_ran_hes.cpp
   src/init_fix_con.cpp
   src/init_fix_like.cpp
   src/init_hes_cross.cpp
   src/init_laplace_obj.cpp
   src/init_laplace_obj_fun.cpp
   src/init_laplace_obj_hes.cpp
   src/fix_con_eval.cpp
   src/fix_con_hes.cpp
   src/fix_con_jac.cpp
   src/fix_like_eval.cpp
   src/fix_like_hes.cpp
   src/fix_like_jac.cpp
   src/logdet_jac.cpp
   src/ran_like_hes.cpp
   src/ran_con_eval.cpp
   src/ran_con_jac.cpp
   src/ran_obj_eval.cpp
   src/ran_obj_jac.cpp
   src/laplace_obj_hes.cpp
   src/update_factor.cpp
}
{xrst_comment --------------------------------------------------------------
   Private Member Variables
}

n_fixed\_
*********
The number of fixed effects is given by
{xrst_spell_off}
{xrst_code cpp} */
   const size_t n_fixed_;
/* {xrst_code}
{xrst_spell_on}
n_random\_
**********
The number of random effects is given by
{xrst_spell_off}
{xrst_code cpp} */
   const size_t n_random_;
/* {xrst_code}
{xrst_spell_on}
quasi_fixed\_
*************
Are we using a quasi-Newton method (or full Newton method)
when :ref:`optimizing fixed effects<optimize_fixed-name>` .
{xrst_spell_off}
{xrst_code cpp} */
   const bool quasi_fixed_;
/* {xrst_code}
{xrst_spell_on}
bool_sparsity\_
***************
If true, use boolean sparsity patterns where possible.
Otherwise, use set sparsity patterns.
{xrst_spell_off}
{xrst_code cpp} */
   const bool bool_sparsity_;
/* {xrst_code}
{xrst_spell_on}

A_rcv\_
*******
contains the random constraint matrix
{xrst_spell_off}
{xrst_code cpp} */
   const d_sparse_rcv A_rcv_;
/* {xrst_code}
{xrst_spell_on}
trace_init\_
************
If true, trace the initialization of cppad_mixed data structures on
standard output. This can be useful for large problems where the initialization
takes a significant amount of time.
{xrst_spell_off}
{xrst_code cpp} */
   const bool trace_init_;
/* {xrst_code}
{xrst_spell_on}

initialize_done\_
*****************
The following flag is false after construction and true after
the corresponding member function is called.
This is the same order as the calls in the file :ref:`initialize-name` :
{xrst_spell_off}
{xrst_code cpp} */
   // only called when n_random_ > 0
   bool                init_ran_con_done_;
   bool                init_ran_like_done_;
   bool                init_ran_jac_done_;
   bool                init_ran_hes_done_;
   bool                init_ldlt_ran_hes_done_;
   bool                init_hes_cross_done_;
   // only called when n_random_ > 0 and quasi_fixed_ is false
   bool                init_laplace_obj_done_;
   bool                init_laplace_obj_fun_done_;
   bool                init_laplace_obj_hes_done_;
   // called in all cases
   bool                init_fix_like_done_;
   bool                init_fix_con_done_;
   // true when all initialization (for this case) is done
   bool                initialize_done_;

/* {xrst_code}
{xrst_spell_on}
cppad_error_handler
*******************
Used to map CppAD error messages to
:ref:`base_class@User Defined Functions@fatal_error` .
{xrst_spell_off}
{xrst_code cpp} */
   CppAD::ErrorHandler cppad_error_handler_;
/* {xrst_code}
{xrst_spell_on}

ran_like_fun\_
**************
If *n_random_*  > 0 and ``init_ran_like_done_`` ,
:ref:`init_ran_like@ran_like_fun_` ,
:ref:`init_ran_like@ran_like_a1fun_` .
are recordings of the user's :ref:`ran_likelihood-name` .
function.
{xrst_spell_off}
{xrst_code cpp} */
   CppAD::ADFun<double>              ran_like_fun_;
   CppAD::ADFun<a1_double, double>   ran_like_a1fun_;
/* {xrst_code}
{xrst_spell_on}
The following objects hold information for computing derivatives
with these ADFun objects:

ran_jac_fun\_
*************
If *n_random_*  > 0 and ``init_ran_jac_done_`` ,
:ref:`init_ran_jac@ran_jac_a1fun_` , and
:ref:`init_ran_jac@ran_jac2hes_rc_` ,
contain the Jacobian of the
:ref:`random likelihood<theory@Random Likelihood, f(theta, u)>`
with respect to the random effects; i.e.
:math:`f_u ( \theta , u )` and the sparsity for
:math:`f_{uu} ( \theta, u )` .
{xrst_spell_off}
{xrst_code cpp} */
   CppAD::ADFun<a1_double, double>  ran_jac_a1fun_;
   sparse_rc                        ran_jac2hes_rc_;
   //
   friend bool ::ran_jac_fun_xam(void);
/* {xrst_code}
{xrst_spell_on}
The row indices in *ran_jac2hes_rc_*
are for just the random effects :math:`u`,
the column indices are for both fixed and random effects
:math:`( \theta , u )`.

ran_hes_fun\_
*************
If *n_random_*  > 0 and ``init_ran_hes_done_`` ,
:ref:`init_ran_hes@ran_hes_uu_rcv_`
contains information for the Hessian of the
:ref:`random likelihood<theory@Random Likelihood, f(theta, u)>`
with respect to the random effects; i.e.
:math:`f_{u,u} ( \theta , u )`.
{xrst_spell_off}
{xrst_code cpp} */
   // sparse Hessian
   d_sparse_rcv                  ran_hes_uu_rcv_;
   //
   // recording of sparse Hessian calculation
   CppAD::ADFun<double>        ran_hes_fun_;
   //
   friend bool ::ran_hes_fun_xam(void);
/* {xrst_code}
{xrst_spell_on}
The row and column indices in *ran_hes_uu_rcv_*
are for just random effects and hence are all less than *n_random* .

ldlt_ran_hes\_
**************
If *n_random_*  > 0 and ``init_ldlt_ran_hes_done_`` ,
``ldlt_ran_hes_`` contains a
:ref:`ldlt_eigen-name` factor for the Hessian of the
:ref:`random likelihood<theory@Random Likelihood, f(theta, u)>`
; i.e.  :math:`f_{u,u} ( \theta , u )`.
{xrst_spell_off}
{xrst_code cpp} */
   CPPAD_MIXED_LDLT_CLASS              ldlt_ran_hes_;
   CppAD::mixed::ldlt_eigen<a1_double> a1_ldlt_ran_hes_;
/* {xrst_code}
{xrst_spell_on}

hes_cross\_
***********
If *n_random_*  > 0 and ``init_hes_cross_done_`` ,
:ref:`init_hes_cross@hes_cross_` contains
information for the cross partials of the Hessian of the
:ref:`random likelihood<theory@Random Likelihood, f(theta, u)>`
; i.e.  :math:`f_{u,\theta} ( \theta , u )`.
{xrst_spell_off}
{xrst_code cpp} */
   CppAD::mixed::sparse_hes_rcv hes_cross_;
   //
   friend bool ::hes_cross_xam(void);
/* {xrst_code}
{xrst_spell_on}

{xrst_comment -------------------------------------------------------------- }

laplace_obj_fun\_
*****************
If *n_random_*  > 0 , quasi_fixed\_ is false, and
``init_laplace_obj_fun_done_`` ,
this is a recording of the second order approximation for the
random part of the Laplace approximation, :math:`H( \beta , \theta , u)`;
see :ref:`init_laplace_obj_fun@laplace_obj_fun_` .
{xrst_spell_off}
{xrst_code cpp} */
   CppAD::ADFun<double>        laplace_obj_fun_;   // for computing H_beta_beta
/* {xrst_code}
{xrst_spell_on}
The following objects hold information for computing derivatives
with this ADFun object:

laplace_obj_hes\_
=================
If *n_random_*  > 0 , quasi_fixed\_ is false, and
``init_laplace_obj_hes_done_`` ,
:ref:`init_laplace_obj_hes@laplace_obj_hes_` contains
information for the Hessian of the
:ref:`Laplace objective<theory@Objective@Laplace Objective, r(theta)>`
{xrst_spell_off}
{xrst_code cpp} */
   CppAD::mixed::sparse_hes_rcv laplace_obj_hes_;
/* {xrst_code}
{xrst_spell_on}
{xrst_comment -------------------------------------------------------------- }

fix_like_fun\_
**************
:ref:`init_fix_like@fix_like_fun_`
is a recording of the fixed part of the likelihood function; see,
:ref:`fix_likelihood-name` .
{xrst_spell_off}
{xrst_code cpp} */
   CppAD::ADFun<double>        fix_like_fun_;     // g(theta)
/* {xrst_code}
{xrst_spell_on}
The following objects hold information for computing derivatives
with this ADFun object:

fix_like_jac\_
==============
:ref:`init_fix_like@fix_like_jac_`
contains information for the Jacobian of the
:ref:`fixed likelihood<theory@Fixed Likelihood, g(theta)>` .
{xrst_spell_off}
{xrst_code cpp} */
   CppAD::mixed::sparse_jac_rcv fix_like_jac_;
/* {xrst_code}
{xrst_spell_on}

fix_like_hes\_
==============
:ref:`init_fix_like@fix_like_hes_`
contains information for the Hessian of the
:ref:`fixed likelihood<theory@Fixed Likelihood, g(theta)>` .
{xrst_spell_off}
{xrst_code cpp} */
   CppAD::mixed::sparse_hes_rcv fix_like_hes_;
/* {xrst_code}
{xrst_spell_on}
{xrst_comment -------------------------------------------------------------- }

fix_con_fun\_
*************
:ref:`init_fix_con@fix_con_fun_`
is a recording of the fixed part of the likelihood function; see,
:ref:`fix_constraint-name` .
{xrst_spell_off}
{xrst_code cpp} */
   CppAD::ADFun<double>        fix_con_fun_;     // c(theta)
/* {xrst_code}
{xrst_spell_on}
The following objects hold information for computing derivatives
with this ADFun object:

fix_con_jac\_
=============
:ref:`init_fix_con@fix_con_jac_`
contains information for the Jacobian of the
constraint function :math:`c ( \theta )`.
{xrst_spell_off}
{xrst_code cpp} */
   CppAD::mixed::sparse_jac_rcv fix_con_jac_;
/* {xrst_code}
{xrst_spell_on}

fix_con_hes\_
=============
If *quasi_fixed* is false,
:ref:`init_fix_con@fix_con_hes_`
contains information for the Hessian of the
:ref:`constraints<fix_constraint-name>` function :math:`c( \theta )`.
The corresponding ADFun object is
:ref:`init_fix_con@fix_con_fun_` .
{xrst_spell_off}
{xrst_code cpp} */
   CppAD::mixed::sparse_hes_rcv fix_con_hes_;
/* {xrst_code}
{xrst_spell_on}
{xrst_comment -------------------------------------------------------------- }
Template Member Functions
*************************
pack
====
See :ref:`pack-name` .
{xrst_spell_off}
{xrst_code cpp} */
   template <class Float_unpack, class Float_pack>
   void pack(
      const CppAD::vector<Float_unpack>& fixed_one  ,
      const CppAD::vector<Float_unpack>& random_vec ,
      CppAD::vector<Float_pack>&         both_vec
   ) const;
   template <class Float_unpack, class Float_pack>
   void pack(
      const CppAD::vector<Float_unpack>& fixed_one  ,
      const CppAD::vector<Float_unpack>& fixed_two  ,
      const CppAD::vector<Float_unpack>& random_vec ,
      CppAD::vector<Float_pack>&         three_vec
   ) const;
/* {xrst_code}
{xrst_spell_on}
unpack
======
See :ref:`unpack-name` .
{xrst_spell_off}
{xrst_code cpp} */
   template <class Float_unpack, class Float_pack>
   void unpack(
      CppAD::vector<Float_unpack>&       fixed_one  ,
      CppAD::vector<Float_unpack>&       random_vec ,
      const CppAD::vector<Float_pack>&   both_vec
   ) const;
   template <class Float_unpack, class Float_pack>
   void unpack(
      CppAD::vector<Float_unpack>&       fixed_one  ,
      CppAD::vector<Float_unpack>&       fixed_two  ,
      CppAD::vector<Float_unpack>&       random_vec ,
      const CppAD::vector<Float_pack>&   three_vec
   ) const;
/* {xrst_code}
{xrst_spell_on}
{xrst_comment -------------------------------------------------------------- }
Initialization Member Functions
*******************************

init_ldlt_ran_hes
=================
See :ref:`init_ldlt_ran_hes-name` .
{xrst_spell_off}
{xrst_code cpp} */
   void init_ldlt_ran_hes(void);
/* {xrst_code}
{xrst_spell_on}

init_fix_con
============
See :ref:`init_fix_con-name` .
{xrst_spell_off}
{xrst_code cpp} */
   void init_fix_con(
      const d_vector& fixed_vec
   );
/* {xrst_code}
{xrst_spell_on}

init_fix_like
=============
See :ref:`init_fix_like-name` .
{xrst_spell_off}
{xrst_code cpp} */
   void init_fix_like(const d_vector& fixed_vec);
/* {xrst_code}
{xrst_spell_on}

init_hes_cross
==============
See :ref:`init_hes_cross-name` .
{xrst_spell_off}
{xrst_code cpp} */
   void init_hes_cross(
      const d_vector& fixed_vec     ,
      const d_vector& random_vec
   );
/* {xrst_code}
{xrst_spell_on}

init_ran_jac
============
See :ref:`init_ran_jac-name` .
{xrst_spell_off}
{xrst_code cpp} */
   void init_ran_jac(
      const d_vector& fixed_vec     ,
      const d_vector& random_vec
   );
/* {xrst_code}
{xrst_spell_on}

init_ran_hes
============
See :ref:`init_ran_hes-name` .
{xrst_spell_off}
{xrst_code cpp} */
   void init_ran_hes(
      const d_vector& fixed_vec     ,
      const d_vector& random_vec
   );
/* {xrst_code}
{xrst_spell_on}

init_laplace_obj
================
See :ref:`init_laplace_obj-name` .
{xrst_spell_off}
{xrst_code cpp} */
   void init_laplace_obj(
      const d_vector& fixed_vec           ,
      const d_vector& random_opt
   );
/* {xrst_code}
{xrst_spell_on}

init_laplace_obj_hes
====================
See :ref:`init_laplace_obj_hes-name` .
{xrst_spell_off}
{xrst_code cpp} */
   void init_laplace_obj_hes(
      const d_vector& fixed_vec     ,
      const d_vector& random_opt
   );
/* {xrst_code}
{xrst_spell_on}

init_ran_like
=============
See :ref:`init_ran_like-name` .
{xrst_spell_off}
{xrst_code cpp} */
   void init_ran_like(
      const d_vector& fixed_vec ,
      const d_vector& random_vec
   );
/* {xrst_code}
{xrst_spell_on}

init_laplace_obj_fun
====================
See :ref:`init_laplace_obj_fun-name` .
{xrst_spell_off}
{xrst_code cpp} */
   void init_laplace_obj_fun(
      const d_vector& fixed_vec ,
      const d_vector& random_opt
   );
/* {xrst_code}
{xrst_spell_on}
{xrst_comment -------------------------------------------------------------- }
Try Member Functions
********************
For each ``try_`` *name* case below,
the public function *name* calls the private function
``try_`` *name* from a ``try`` block
with a corresponding ``catch`` that maps a
``cppad_mixed`` :ref:`exception-name` to a
:ref:`base_class@User Defined Functions@fatal_error` call.

try_initialize
==============
Called by public :ref:`base_class@initialize`
{xrst_spell_off}
{xrst_code cpp} */
   std::map<std::string, size_t> try_initialize(
      const d_vector&  fixed_vec  ,
      const d_vector&  random_vec
   );
/* {xrst_code}
{xrst_spell_on}
try_optimize_random
===================
Called by public :ref:`base_class@optimize_random`
{xrst_spell_off}
{xrst_code cpp} */
   d_vector try_optimize_random(
      const std::string& options      ,
      const d_vector&    fixed_vec    ,
      const d_vector&    random_lower ,
      const d_vector&    random_upper ,
      const d_vector&    random_in
   );
/* {xrst_code}
{xrst_spell_on}
try_optimize_fixed
==================
Called by public :ref:`base_class@optimize_fixed`
{xrst_spell_off}
{xrst_code cpp} */
   CppAD::mixed::fixed_solution try_optimize_fixed(
      const std::string& fixed_ipopt_options   ,
      const std::string& random_ipopt_options  ,
      const d_vector&    fixed_lower           ,
      const d_vector&    fixed_upper           ,
      const d_vector&    fix_constraint_lower  ,
      const d_vector&    fix_constraint_upper  ,
      const d_vector&    fixed_scale           ,
      const d_vector&    fixed_in              ,
      const d_vector&    random_lower          ,
      const d_vector&    random_upper          ,
      const d_vector&    random_in             ,
      const CppAD::mixed::warm_start_struct& warm_start
   );
/* {xrst_code}
{xrst_spell_on}
try_hes_fixed_obj
=================
Called by public
:ref:`base_class@hes_fixed_obj` and
:ref:`information_mat<base_class@information_mat, Deprecated 2020-03-22>` .
{xrst_spell_off}
{xrst_code cpp} */
   d_sparse_rcv try_hes_fixed_obj(
      const d_vector& fixed_vec      ,
      const d_vector& random_opt
   );
/* {xrst_code}
{xrst_spell_on}
try_hes_random_obj
==================
Called by public :ref:`base_class@hes_random_obj`
{xrst_spell_off}
{xrst_code cpp} */
   d_sparse_rcv try_hes_random_obj(
      const d_vector& fixed_vec      ,
      const d_vector& random_vec
   );
/* {xrst_code}
{xrst_spell_on}
try_sample_fixed
================
Called by public :ref:`base_class@sample_fixed`
{xrst_spell_off}
{xrst_code cpp} */
   std::string try_sample_fixed(
      d_vector&                            sample               ,
      const d_sparse_rcv&                  information_rcv      ,
      const CppAD::mixed::fixed_solution&  solution             ,
      const d_vector&                      fixed_lower          ,
      const d_vector&                      fixed_upper
   );
/* {xrst_code}
{xrst_spell_on}
sample_random
=============
Called by public :ref:`base_class@sample_random`
{xrst_spell_off}
{xrst_code cpp} */
   std::string try_sample_random(
      d_vector&             sample               ,
      const std::string&    random_ipopt_options ,
      const d_vector&       fixed_vec            ,
      const d_vector&       random_lower         ,
      const d_vector&       random_upper         ,
      const d_vector&       random_in
   );
/* {xrst_code}
{xrst_spell_on}
{xrst_comment -------------------------------------------------------------- }
Other Member Functions
**********************

fix_con_eval
============
See :ref:`fix_con_eval-name`
{xrst_spell_off}
{xrst_code cpp} */
   d_vector fix_con_eval(const d_vector& fixed_vec);
   friend bool ::fix_con_eval_xam(void);
/* {xrst_code}
{xrst_spell_on}

fix_con_jac
===========
See :ref:`fix_con_jac-name`
{xrst_spell_off}
{xrst_code cpp} */
   void fix_con_jac(
      const d_vector&        fixed_vec   ,
      CppAD::vector<size_t>& row_out     ,
      CppAD::vector<size_t>& col_out     ,
      d_vector&              val_out
   );
   friend bool ::fix_con_jac_xam(void);
/* {xrst_code}
{xrst_spell_on}

fix_con_hes
===========
See :ref:`fix_con_hes-name`
{xrst_spell_off}
{xrst_code cpp} */
   void fix_con_hes(
      const d_vector&        fixed_vec   ,
      const d_vector&        weight      ,
      CppAD::vector<size_t>& row_out     ,
      CppAD::vector<size_t>& col_out     ,
      d_vector&              val_out
   );
   friend bool ::fix_con_hes_xam(void);
/* {xrst_code}
{xrst_spell_on}

fix_like_eval
=============
See :ref:`fix_like_eval-name`
{xrst_spell_off}
{xrst_code cpp} */
   d_vector fix_like_eval(const d_vector& fixed_vec);
   friend bool ::fix_like_eval_xam(void);
/* {xrst_code}
{xrst_spell_on}

fix_like_jac
============
See :ref:`fix_like_jac-name`
{xrst_spell_off}
{xrst_code cpp} */
   void fix_like_jac(
      const d_vector&        fixed_vec   ,
      CppAD::vector<size_t>& row_out     ,
      CppAD::vector<size_t>& col_out     ,
      d_vector&              val_out
   );
   friend bool ::fix_like_jac_xam(void);
/* {xrst_code}
{xrst_spell_on}

fix_like_hes
============
See :ref:`fix_like_hes-name`
{xrst_spell_off}
{xrst_code cpp} */
   void fix_like_hes(
      const d_vector&        fixed_vec   ,
      const d_vector&        weight      ,
      CppAD::vector<size_t>& row_out     ,
      CppAD::vector<size_t>& col_out     ,
      d_vector&              val_out
   );
   friend bool ::fix_like_hes_xam(void);
/* {xrst_code}
{xrst_spell_on}

logdet_jac
==========
See :ref:`logdet_jac-name`
{xrst_spell_off}
{xrst_code cpp} */
   void logdet_jac(
      const d_vector& fixed_vec  ,
      const d_vector& random_vec ,
      d_vector&       logdet_fix ,
      d_vector&       logdet_ran
   );
   friend bool ::logdet_jac_xam(void);
/* {xrst_code}
{xrst_spell_on}

ran_con_eval
============
See :ref:`ran_con_eval-name`
{xrst_spell_off}
{xrst_code cpp} */
   void ran_con_eval(
      const d_vector& random_vec ,
      d_vector&       Au
   );
   friend bool ::ran_con_eval_xam(void);
/* {xrst_code}
{xrst_spell_on}

ran_con_jac
===========
See :ref:`ran_con_jac-name`
{xrst_spell_off}
{xrst_code cpp} */
   void ran_con_jac(
      const d_vector&                fixed_vec  ,
      const d_vector&                random_vec ,
      d_sparse_rcv&                  jac_info
   );
   friend bool ::ran_con_jac_xam(void);
   friend bool ::sample_fixed_1(void);
/* {xrst_code}
{xrst_spell_on}

ran_obj_eval
============
See :ref:`ran_obj_eval-name`
{xrst_spell_off}
{xrst_code cpp} */
   double ran_obj_eval(
      const d_vector& fixed_vec  ,
      const d_vector& random_vec
   );
   friend bool ::ran_obj_eval_xam(void);
   friend bool ::laplace_obj_fun(void);
   friend bool ::laplace_obj_tst(void);
/* {xrst_code}
{xrst_spell_on}

ran_obj_jac
===========
See :ref:`ran_obj_jac-name`
{xrst_spell_off}
{xrst_code cpp} */
   void ran_obj_jac(
      const d_vector& fixed_vec  ,
      const d_vector& random_vec ,
      d_vector&       r_fixed
   );
   friend bool ::ran_obj_jac_xam(void);
   friend bool ::der_var_hes(void);
   friend bool ::delta_ran_obj(void);
/* {xrst_code}
{xrst_spell_on}

laplace_obj_hes
===============
See :ref:`laplace_obj_hes-name`
{xrst_spell_off}
{xrst_code cpp} */
   void laplace_obj_hes(
      const d_vector&         fixed_vec   ,
      const d_vector&         random_vec  ,
      const d_vector&         weight      ,
      CppAD::vector<size_t>&  row_out     ,
      CppAD::vector<size_t>&  col_out     ,
      d_vector&               val_out
   );
   friend bool ::laplace_obj_hes_xam(void);
   friend bool ::laplace_obj_hes(void);
/* {xrst_code}
{xrst_spell_on}

update_factor
=============
See :ref:`update_factor-name`
{xrst_spell_off}
{xrst_code cpp} */
   void update_factor(
      const d_vector&         fixed_vec   ,
      const d_vector&         random_vec
   );
   friend bool ::update_factor_xam(void);
/* {xrst_code}
{xrst_spell_on}

{xrst_end private_base_class}
-------------------------------------------------------------------------------
*/
# endif
