// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-25 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_BASE_CLASS_HPP
# define CPPAD_MIXED_BASE_CLASS_HPP

// cppad_mixed.hpp: public member fucntions
public:
/*
{xrst_begin base_class}

cppad_mixed: Public Declarations
################################
These ``cppad_mixed`` class declarations are ``public`` .
They are part of the user API and can be used by a derived class object
:ref:`derived_ctor@mixed_object` .

Cppad Mixed Types
*****************
The following Cppad Mixed types are extended into the
``cppad_mixed`` class:
{xrst_spell_off}
{xrst_code cpp} */
   // scalar types
   typedef CppAD::mixed::a1_double      a1_double;
   // vector types
   typedef CppAD::mixed::s_vector       s_vector;
   typedef CppAD::mixed::d_vector       d_vector;
   typedef CppAD::mixed::a1_vector      a1_vector;
   // sparse types
   typedef CppAD::mixed::sparse_rc      sparse_rc;
   typedef CppAD::mixed::d_sparse_rcv   d_sparse_rcv;
   typedef CppAD::mixed::a1_sparse_rcv  a1_sparse_rcv;
/* {xrst_code}
{xrst_spell_on}

User Defined Functions
**********************
The following are ``cppad_mixed`` pure virtual functions.
Each one has a default definition that may be replaced by
the user's derived class:

ran_likelihood
==============
This function is necessary if there are random effects in the model.
{xrst_spell_off}
{xrst_code cpp} */
   virtual a1_vector ran_likelihood(
      const a1_vector& fixed_vec  ,
      const a1_vector& random_vec )
   {  return a1_vector(0); }
/* {xrst_code}
{xrst_spell_on}
See :ref:`ran_likelihood-name` .

fix_likelihood
==============
This function should be used if there is a prior on the fixed effects,
or there is data that does not depend on the random effects.
{xrst_spell_off}
{xrst_code cpp} */
   virtual a1_vector fix_likelihood(
      const a1_vector& fixed_vec )
   {  return a1_vector(0); }
/* {xrst_code}
{xrst_spell_on}
See :ref:`fix_likelihood-name` .

fix_constraint
==============
This function is used to define constraints
that only depend on the fixed effects.
{xrst_spell_off}
{xrst_code cpp} */
   virtual a1_vector fix_constraint(
      const a1_vector& fixed_vec )
   {  return a1_vector(0); }
/* {xrst_code}
{xrst_spell_on}
See :ref:`fix_constraint-name` .

fatal_error
===========
This routine displays an error message and then exits the program.
Its default definition below can be replaced by a user definition.
Note that if ``NDEBUG`` is not defined, this generates an assert,
otherwise it exits.
{xrst_spell_off}
{xrst_code cpp} */
   virtual void fatal_error(const std::string& error_message)
   {  std::cerr << "cppad_mixed error: " << error_message << std::endl;
      assert(false);
      exit(1);
   }

/* {xrst_code}
{xrst_spell_on}

warning
=======
This routine displays a warning message and then returns.
Its default definition below can be replaced by a user definition.
{xrst_spell_off}
{xrst_code cpp} */
   virtual void warning(const std::string& warning_message)
   {  std::cerr << "cppad_mixed warning: " << warning_message << std::endl;
   }
/* {xrst_code}
{xrst_spell_on}

constructor
***********
:ref:`derived_ctor-title`.
{xrst_spell_off}
{xrst_code cpp} */
   // BEGIN_CPPAD_MIXED_CTOR
   cppad_mixed(
      size_t               n_fixed       ,
      size_t               n_random      ,
      bool                 quasi_fixed   = false           ,
      bool                 bool_sparsity = false           ,
      const d_sparse_rcv&  A_rcv         = d_sparse_rcv()  ,
      bool                 trace_init    = false
   );
   // END_CPPAD_MIXED_CTOR
/* {xrst_code}
{xrst_spell_on}

destructor
**********
{xrst_spell_off}
{xrst_code cpp} */
   virtual ~cppad_mixed(void);
/* {xrst_code}
{xrst_spell_on}

initialize
**********
:ref:`initialize-title`
{xrst_spell_off}
{xrst_code cpp} */
   std::map<std::string, size_t> initialize(
      const d_vector&  fixed_vec   ,
      const d_vector&  random_vec
   );
/* {xrst_code}
{xrst_spell_on}
optimize_random
***************
:ref:`optimize_random-title`.
{xrst_spell_off}
{xrst_code cpp} */
   d_vector optimize_random(
      const std::string& options      ,
      const d_vector&    fixed_vec    ,
      const d_vector&    random_lower ,
      const d_vector&    random_upper ,
      const d_vector&    random_in
   );
/* {xrst_code}
{xrst_spell_on}
optimize_fixed
**************
:ref:`optimize_fixed-title`.
{xrst_spell_off}
{xrst_code cpp} */
   CppAD::mixed::fixed_solution optimize_fixed(
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
      const CppAD::mixed::warm_start_struct&  warm_start =
         CppAD::mixed::warm_start_struct()
   );
/* {xrst_code}
{xrst_spell_on}
hes_fixed_obj
*************
:ref:`hes_fixed_obj-title`.
{xrst_spell_off}
{xrst_code cpp} */
   d_sparse_rcv hes_fixed_obj(
      const d_vector& fixed_vec  ,
      const d_vector& random_opt
   );
/* {xrst_code}
{xrst_spell_on}
hes_random_obj
**************
:ref:`hes_random_obj-title`.
{xrst_spell_off}
{xrst_code cpp} */
   d_sparse_rcv hes_random_obj(
      const d_vector& fixed_vec  ,
      const d_vector& random_vec
   );
/* {xrst_code}
{xrst_spell_on}
sample_fixed
************
:ref:`sample_fixed-title`.
{xrst_spell_off}
{xrst_code cpp} */
   std::string sample_fixed(
      d_vector&                            sample               ,
      const d_sparse_rcv&                  information_rcv      ,
      const CppAD::mixed::fixed_solution&  solution             ,
      const d_vector&                      fixed_lower          ,
      const d_vector&                      fixed_upper          ,
      double&                              rcond
   );
   std::string sample_fixed(
      d_vector&                            sample               ,
      const d_sparse_rcv&                  information_rcv      ,
      const CppAD::mixed::fixed_solution&  solution             ,
      const d_vector&                      fixed_lower          ,
      const d_vector&                      fixed_upper          )
   {  double rcond;
      return sample_fixed(
         sample, information_rcv, solution, fixed_lower, fixed_upper, rcond
      );
   }
/* {xrst_code}
{xrst_spell_on}
sample_random
*************
:ref:`sample_random-title`.
{xrst_spell_off}
{xrst_code cpp} */
   std::string sample_random(
      d_vector&            sample               ,
      const std::string&   random_ipopt_options ,
      const d_vector&      fixed_vec            ,
      const d_vector&      random_lower         ,
      const d_vector&      random_upper         ,
      const d_vector&      random_in
   );
/* {xrst_code}
{xrst_spell_on}
information_mat, Deprecated 2020-03-22
**************************************
:ref:`information_mat-title`.
{xrst_spell_off}
{xrst_code cpp} */
   d_sparse_rcv information_mat(
      const CppAD::mixed::fixed_solution&  solution             ,
      const d_vector&                      random_opt
   );
/* {xrst_code}
{xrst_spell_on}

Contents
********
{xrst_toc_table
   src/derived_ctor.cpp
   src/ran_likelihood.xrst
   src/fix_likelihood.xrst
   src/fix_constraint.xrst
   src/initialize.cpp
   src/optimize_random.cpp
   src/optimize_fixed.cpp
   src/eigen/hes_fixed_obj.cpp
   src/eigen/hes_random_obj.cpp
   src/sample_fixed.cpp
   src/sample_random.cpp
   src/eigen/information_mat.cpp
}

{xrst_end base_class}
*/

# endif
