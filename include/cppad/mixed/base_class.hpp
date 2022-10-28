// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_BASE_CLASS_HPP
# define CPPAD_MIXED_BASE_CLASS_HPP

// cppad_mixed.hpp: public member fucntions
public:
/*
$begin base_class$$
$spell
   rcv
   jac
   ldlt_eigen
   std
   cerr
   endl
   CppAD
   init
   ran_obj
   cppad
   obj
   const
   vec
   typedef
   CppAD
   ctor
   std
   bool
   hes
$$

$section cppad_mixed: Public Declarations$$
These $code cppad_mixed$$ class declarations are $code public$$.
They are part of the user API and can be used by a derived class object
$cref/mixed_object/derived_ctor/mixed_object/$$.

$head Cppad Mixed Types$$
The following Cppad Mixed types are extended into the
$code cppad_mixed$$ class:
$srccode%cpp% */
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
/* %$$

$head User Defined Functions$$
The following are $code cppad_mixed$$ pure virtual functions.
Each one has a default definition that may be replaced by
the user's derived class:

$subhead ran_likelihood$$
This function is necessary if there are random effects in the model.
$srccode%cpp% */
   virtual a1_vector ran_likelihood(
      const a1_vector& fixed_vec  ,
      const a1_vector& random_vec )
   {  return a1_vector(0); }
/* %$$
See $cref ran_likelihood$$.

$subhead fix_likelihood$$
This function should be used if there is a prior on the fixed effects,
or there is data that does not depend on the random effects.
$srccode%cpp% */
   virtual a1_vector fix_likelihood(
      const a1_vector& fixed_vec )
   {  return a1_vector(0); }
/* %$$
See $cref fix_likelihood$$.

$subhead fix_constraint$$
This function is used to define constraints
that only depend on the fixed effects.
$srccode%cpp% */
   virtual a1_vector fix_constraint(
      const a1_vector& fixed_vec )
   {  return a1_vector(0); }
/* %$$
See $cref fix_constraint$$.

$subhead fatal_error$$
This routine displays an error message and then exits the program.
Its default definition below can be replaced by a user definition.
Note that if $code NDEBUG$$ is not defined, this generates an assert,
otherwise it exits.
$srccode%cpp% */
   virtual void fatal_error(const std::string& error_message)
   {  std::cerr << "cppad_mixed error: " << error_message << std::endl;
      assert(false);
      exit(1);
   }

/* %$$

$subhead warning$$
This routine displays a warning message and then returns.
Its default definition below can be replaced by a user definition.
$srccode%cpp% */
   virtual void warning(const std::string& warning_message)
   {  std::cerr << "cppad_mixed warning: " << warning_message << std::endl;
   }
/* %$$

$head constructor$$
$cref derived_ctor$$, $title derived_ctor$$.
$srccode%cpp% */
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
/* %$$

$head destructor$$
$srccode%cpp% */
   virtual ~cppad_mixed(void);
/* %$$

$head initialize$$
$cref initialize$$, $title initialize$$
$srccode%cpp% */
   std::map<std::string, size_t> initialize(
      const d_vector&  fixed_vec   ,
      const d_vector&  random_vec
   );
/* %$$
$head optimize_random$$
$cref optimize_random$$, $title optimize_random$$.
$srccode%cpp% */
   d_vector optimize_random(
      const std::string& options      ,
      const d_vector&    fixed_vec    ,
      const d_vector&    random_lower ,
      const d_vector&    random_upper ,
      const d_vector&    random_in
   );
/* %$$
$head optimize_fixed$$
$cref optimize_fixed$$, $title optimize_fixed$$.
$srccode%cpp% */
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
/* %$$
$head hes_fixed_obj$$
$cref hes_fixed_obj$$, $title hes_fixed_obj$$.
$srccode%cpp% */
   d_sparse_rcv hes_fixed_obj(
      const d_vector& fixed_vec  ,
      const d_vector& random_opt
   );
/* %$$
$head hes_random_obj$$
$cref hes_random_obj$$, $title hes_random_obj$$.
$srccode%cpp% */
   d_sparse_rcv hes_random_obj(
      const d_vector& fixed_vec  ,
      const d_vector& random_vec
   );
/* %$$
$head sample_fixed$$
$cref sample_fixed$$, $title sample_fixed$$.
$srccode%cpp% */
   std::string sample_fixed(
      d_vector&                            sample               ,
      const d_sparse_rcv&                  information_rcv      ,
      const CppAD::mixed::fixed_solution&  solution             ,
      const d_vector&                      fixed_lower          ,
      const d_vector&                      fixed_upper
   );
/* %$$
$head sample_random$$
$cref sample_random$$, $title sample_random$$.
$srccode%cpp% */
   std::string sample_random(
      d_vector&            sample               ,
      const std::string&   random_ipopt_options ,
      const d_vector&      fixed_vec            ,
      const d_vector&      random_lower         ,
      const d_vector&      random_upper         ,
      const d_vector&      random_in
   );
/* %$$
$head information_mat, Deprecated 2020-03-22$$
$cref information_mat$$, $title information_mat$$.
$srccode%cpp% */
   d_sparse_rcv information_mat(
      const CppAD::mixed::fixed_solution&  solution             ,
      const d_vector&                      random_opt
   );
/* %$$
$childtable%src/derived_ctor.cpp
   %src/ran_likelihood.omh
   %src/fix_likelihood.omh
   %src/fix_constraint.omh
   %src/initialize.cpp
   %src/optimize_random.cpp
   %src/optimize_fixed.cpp
   %src/eigen/hes_fixed_obj.cpp
   %src/eigen/hes_random_obj.cpp
   %src/sample_fixed.cpp
   %src/sample_random.cpp
   %src/eigen/information_mat.cpp
%$$


$end
*/

# endif
