/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-20 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_PRIVATE_BASE_CLASS_HPP
# define CPPAD_MIXED_PRIVATE_BASE_CLASS_HPP

// cppad_mixed.hpp: private member functions and variables
private:
/*
------------------------------------------------------------------------------
$begin private_base_class$$
$spell
	rc
	rcv
	ldlt_eigen
	objcon
	eigen
	chol
	Cholesky
	CppAD
	init
	ran_obj
	var
	cppad
	hes hes
	obj
	jac
	Jacobians
	jacobian
	hes
	eval
	typedef
	CppAD
	vec
	const
	bool
	xam
	uu
	logdet
	Au
$$

$section cppad_mixed: Private Declarations$$
These $code cppad_mixed$$ class declarations are $code private$$.
They are $bold not$$ part of the user API,
they may change with time, and they
can $bold not$$ be used by a derived class object
$cref/mixed_object/derived_ctor/mixed_object/$$.

$childtable%include/cppad/mixed/pack.hpp
	%include/cppad/mixed/unpack.hpp

	%src/init_ran_like.cpp
	%src/init_ran_jac.cpp
	%src/init_ran_hes.cpp
	%src/init_ldlt_ran_hes.cpp
	%src/init_fix_con.cpp
	%src/init_fix_like.cpp
	%src/init_hes_cross.cpp
	%src/init_laplace_obj.cpp
	%src/init_laplace_obj_fun.cpp
	%src/init_laplace_obj_hes.cpp

	%src/fix_con_eval.cpp
	%src/fix_con_hes.cpp
	%src/fix_con_jac.cpp
	%src/fix_like_eval.cpp
	%src/fix_like_hes.cpp
	%src/fix_like_jac.cpp
	%src/logdet_jac.cpp
	%src/ran_like_hes.cpp
	%src/ran_con_eval.cpp
	%src/ran_con_jac.cpp
	%src/ran_obj_eval.cpp
	%src/ran_obj_jac.cpp
	%src/laplace_obj_hes.cpp
	%src/update_factor.cpp
%$$
$comment -------------------------------------------------------------------
	Private Member Variables
$$

$head n_fixed_$$
The number of fixed effects is given by
$srccode%cpp% */
	const size_t n_fixed_;
/* %$$
$head n_random_$$
The number of random effects is given by
$srccode%cpp% */
	const size_t n_random_;
/* %$$
$head quasi_fixed_$$
Are we using a quasi-Newton method (or full Newton method)
when $cref/optimizing fixed effects/optimize_fixed/$$.
$srccode%cpp% */
	const bool quasi_fixed_;
/* %$$
$head bool_sparsity_$$
f true, use boolean sparsity patterns where possible.
Otherwise, use set sparsity patterns.
$srccode%cpp% */
	const bool bool_sparsity_;
/* %$$

$head A_rcv_$$
contains the random constraint matrix
$srccode%cpp% */
	const d_sparse_rcv A_rcv_;
/* %$$

$head initialize_done_$$
The following flag is false after construction and true after
the corresponding member function is called.
This is the same order as the calls in the file $cref initialize$$:
$srccode%cpp% */
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

/* %$$
$head cppad_error_handler$$
Used to map CppAD error messages to
$cref/fatal_error/base_class/User Defined Functions/fatal_error/$$.
$srccode%cpp% */
	CppAD::ErrorHandler cppad_error_handler_;
/* %$$

$head ran_like_fun_$$
If $icode%n_random_% > 0%$$ and $code init_ran_like_done_$$,
$cref/ran_like_fun_/init_ran_like/ran_like_fun_/$$,
$cref/ran_like_a1fun_/init_ran_like/ran_like_a1fun_/$$.
are recordings of the user's $cref ran_likelihood$$.
function.
$srccode%cpp% */
	CppAD::ADFun<double>              ran_like_fun_;
	CppAD::ADFun<a1_double, double>   ran_like_a1fun_;
/* %$$
The following objects hold information for computing derivatives
with these ADFun objects:

$head ran_jac_fun_$$
If $icode%n_random_% > 0%$$ and $code init_ran_jac_done_$$,
$cref/ran_jac_a1fun_/init_ran_jac/ran_jac_a1fun_/$$, and
$cref/ran_jac2hes_rc_/init_ran_jac/ran_jac2hes_rc_/$$,
contain the Jacobian of the
$cref/random likelihood
	/theory
	/Random Likelihood, f(theta, u)
/$$
with respect to the random effects; i.e.
$latex f_u ( \theta , u )$$ and the sparsity for
$latex f_{uu} ( \theta, u )$$ .
$srccode%cpp% */
	CppAD::ADFun<a1_double, double>  ran_jac_a1fun_;
	sparse_rc                        ran_jac2hes_rc_;
	//
	friend bool ::ran_jac_fun_xam(void);
/* %$$
The row indices in $icode ran_jac2hes_rc_$$
are for just the random effects $latex u$$,
the column indices are for both fixed and random effects
$latex ( \theta , u )$$.

$head ran_hes_fun_$$
If $icode%n_random_% > 0%$$ and $code init_ran_hes_done_$$,
$cref/ran_hes_uu_rcv_/init_ran_hes/ran_hes_uu_rcv_/$$
contains information for the Hessian of the
$cref/random likelihood
	/theory
	/Random Likelihood, f(theta, u)
/$$
with respect to the random effects; i.e.
$latex f_{u,u} ( \theta , u )$$.
$srccode%cpp% */
	// sparse Hessian
	d_sparse_rcv                  ran_hes_uu_rcv_;
	//
	// recording of sparse Hessian calculation
	CppAD::ADFun<double>        ran_hes_fun_;
	//
	friend bool ::ran_hes_fun_xam(void);
/* %$$
The row and column indices in $icode ran_hes_uu_rcv_$$
are for just random effects and hence are all less than $icode n_random$$.

$head ldlt_ran_hes_$$
If $icode%n_random_% > 0%$$ and $code init_ldlt_ran_hes_done_$$,
$code ldlt_ran_hes_$$ contains a
$cref ldlt_eigen$$ factor for the Hessian of the
$cref/random likelihood
	/theory
	/Random Likelihood, f(theta, u)
/$$
; i.e.  $latex f_{u,u} ( \theta , u )$$.
$srccode%cpp% */
	CPPAD_MIXED_LDLT_CLASS              ldlt_ran_hes_;
	CppAD::mixed::ldlt_eigen<a1_double> a1_ldlt_ran_hes_;
/* %$$

$head hes_cross_$$
If $icode%n_random_% > 0%$$ and $code init_hes_cross_done_$$,
$cref/hes_cross_/init_hes_cross/hes_cross_/$$ contains
information for the cross partials of the Hessian of the
$cref/random likelihood
	/theory
	/Random Likelihood, f(theta, u)
/$$
; i.e.  $latex f_{u,\theta} ( \theta , u )$$.
$srccode%cpp% */
	CppAD::mixed::sparse_hes_rcv hes_cross_;
	//
	friend bool ::hes_cross_xam(void);
/* %$$

$comment ------------------------------------------------------------------- $$

$head laplace_obj_fun_$$
If $icode%n_random_% > 0%$$, quasi_fixed_ is false, and
$code init_laplace_obj_fun_done_$$,
this is a recording of the second order approximation for the
random part of the Laplace approximation, $latex H( \beta , \theta , u)$$;
see $cref/laplace_obj_fun_/init_laplace_obj_fun/laplace_obj_fun_/$$.
$srccode%cpp% */
	CppAD::ADFun<double>        laplace_obj_fun_;   // for computing H_beta_beta
/* %$$
The following objects hold information for computing derivatives
with this ADFun object:

$subhead laplace_obj_hes_$$
If $icode%n_random_% > 0%$$, quasi_fixed_ is false, and
$code init_laplace_obj_hes_done_$$,
$cref/laplace_obj_hes_/init_laplace_obj_hes/laplace_obj_hes_/$$ contains
information for the Hessian of the
$cref/Laplace objective
	/theory
	/Objective
	/Laplace Objective, r(theta)
/$$
$srccode%cpp% */
	CppAD::mixed::sparse_hes_rcv laplace_obj_hes_;
/* %$$
$comment ------------------------------------------------------------------- $$

$head fix_like_fun_$$
$cref/fix_like_fun_/init_fix_like/fix_like_fun_/$$
is a recording of the fixed part of the likelihood function; see,
$cref fix_likelihood$$.
$srccode%cpp% */
	CppAD::ADFun<double>        fix_like_fun_;     // g(theta)
/* %$$
The following objects hold information for computing derivatives
with this ADFun object:

$subhead fix_like_jac_$$
$cref/fix_like_jac_/init_fix_like/fix_like_jac_/$$
contains information for the Jacobian of the
$cref/fixed likelihood/theory/Fixed Likelihood, g(theta)/$$.
$srccode%cpp% */
	CppAD::mixed::sparse_jac_rcv fix_like_jac_;
/* %$$

$subhead fix_like_hes_$$
$cref/fix_like_hes_/init_fix_like/fix_like_hes_/$$
contains information for the Hessian of the
$cref/fixed likelihood/theory/Fixed Likelihood, g(theta)/$$.
$srccode%cpp% */
	CppAD::mixed::sparse_hes_rcv fix_like_hes_;
/* %$$
$comment ------------------------------------------------------------------- $$

$head fix_con_fun_$$
$cref/fix_con_fun_/init_fix_con/fix_con_fun_/$$
is a recording of the fixed part of the likelihood function; see,
$cref fix_constraint$$.
$srccode%cpp% */
	CppAD::ADFun<double>        fix_con_fun_;     // c(theta)
/* %$$
The following objects hold information for computing derivatives
with this ADFun object:

$subhead fix_con_jac_$$
$cref/fix_con_jac_/init_fix_con/fix_con_jac_/$$
contains information for the Jacobian of the
constraint function $latex c ( \theta )$$.
$srccode%cpp% */
	CppAD::mixed::sparse_jac_rcv fix_con_jac_;
/* %$$

$subhead fix_con_hes_$$
If $icode quasi_fixed$$ is false,
$cref/fix_con_hes_/init_fix_con/fix_con_hes_/$$
contains information for the Hessian of the
$cref/constraints/fix_constraint/$$ function $latex c( \theta )$$.
The corresponding ADFun object is
$cref/fix_con_fun_/init_fix_con/fix_con_fun_/$$.
$srccode%cpp% */
	CppAD::mixed::sparse_hes_rcv fix_con_hes_;
/* %$$
$comment ------------------------------------------------------------------- $$
$head Template Member Functions$$
$subhead pack$$
See $cref pack$$.
$srccode%cpp% */
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
/* %$$
$subhead unpack$$
See $cref unpack$$.
$srccode%cpp% */
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
/* %$$
$comment ------------------------------------------------------------------- $$
$head Initialization Member Functions$$

$subhead init_ldlt_ran_hes$$
See $cref init_ldlt_ran_hes$$.
$srccode%cpp% */
	void init_ldlt_ran_hes(void);
/* %$$

$subhead init_fix_con$$
See $cref init_fix_con$$.
$srccode%cpp% */
	void init_fix_con(
		const d_vector& fixed_vec
	);
/* %$$

$subhead init_fix_like$$
See $cref init_fix_like$$.
$srccode%cpp% */
	void init_fix_like(const d_vector& fixed_vec);
/* %$$

$subhead init_hes_cross$$
See $cref init_hes_cross$$.
$srccode%cpp% */
	void init_hes_cross(
		const d_vector& fixed_vec     ,
		const d_vector& random_vec
	);
/* %$$

$subhead init_ran_jac$$
See $cref init_ran_jac$$.
$srccode%cpp% */
	void init_ran_jac(
		const d_vector& fixed_vec     ,
		const d_vector& random_vec
	);
/* %$$

$subhead init_ran_hes$$
See $cref init_ran_hes$$.
$srccode%cpp% */
	void init_ran_hes(
		const d_vector& fixed_vec     ,
		const d_vector& random_vec
	);
/* %$$

$subhead init_laplace_obj$$
See $cref init_laplace_obj$$.
$srccode%cpp% */
	void init_laplace_obj(
		const d_vector& fixed_vec           ,
		const d_vector& random_in           ,
		const d_vector& random_lower        ,
		const d_vector& random_upper        ,
		const std::string& random_options
	);
/* %$$

$subhead init_laplace_obj_hes$$
See $cref init_laplace_obj_hes$$.
$srccode%cpp% */
	void init_laplace_obj_hes(
		const d_vector& fixed_vec     ,
		const d_vector& random_vec
	);
/* %$$

$subhead init_ran_like$$
See $cref init_ran_like$$.
$srccode%cpp% */
	void init_ran_like(
		const d_vector& fixed_vec ,
		const d_vector& random_vec
	);
/* %$$

$subhead init_laplace_obj_fun$$
See $cref init_laplace_obj_fun$$.
$srccode%cpp% */
	void init_laplace_obj_fun(
		const d_vector& fixed_vec ,
		const d_vector& random_vec
	);
/* %$$
$comment ------------------------------------------------------------------- $$
$head Try Member Functions$$
For each $codei%try_%name%$$ case below,
the public function $icode name$$ calls the private function
$codei%try_%name%$$ from a $code try$$ block
with a corresponding $code catch$$ that maps a
$code cppad_mixed$$ $cref exception$$ to a
$cref/fatal_error/base_class/User Defined Functions/fatal_error/$$ call.

$subhead try_initialize$$
Called by public $cref/initialize/base_class/initialize/$$
$srccode%cpp% */
	std::map<std::string, size_t> try_initialize(
		const d_vector&  fixed_vec  ,
		const d_vector&  random_vec
	);
/* %$$
$subhead try_optimize_random$$
Called by public $cref/optimize_random/base_class/optimize_random/$$
$srccode%cpp% */
	d_vector try_optimize_random(
		const std::string& options      ,
		const d_vector&    fixed_vec    ,
		const d_vector&    random_lower ,
		const d_vector&    random_upper ,
		const d_vector&    random_in
	);
/* %$$
$subhead try_optimize_fixed$$
Called by public $cref/optimize_fixed/base_class/optimize_fixed/$$
$srccode%cpp% */
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
		const d_vector&    random_in
	);
/* %$$
$subhead try_hes_fixed_obj$$
Called by public
$cref/hes_fixed_obj/base_class/hes_fixed_obj/$$ and
$cref/information_mat/base_class/information_mat, Deprecated 2020-03-22/$$.
$srccode%cpp% */
	d_sparse_rcv try_hes_fixed_obj(
		const d_vector& fixed_vec      ,
		const d_vector& random_opt
	);
/* %$$
$subhead try_hes_random_obj$$
Called by public $cref/hes_random_obj/base_class/hes_random_obj/$$
$srccode%cpp% */
	d_sparse_rcv try_hes_random_obj(
		const d_vector& fixed_vec      ,
		const d_vector& random_vec
	);
/* %$$
$subhead try_sample_fixed$$
Called by public $cref/sample_fixed/base_class/sample_fixed/$$
$srccode%cpp% */
	std::string try_sample_fixed(
		d_vector&                            sample               ,
		const d_sparse_rcv&                  information_rcv      ,
		const CppAD::mixed::fixed_solution&  solution             ,
		const d_vector&                      fixed_lower          ,
		const d_vector&                      fixed_upper
	);
/* %$$
$subhead sample_random$$
Called by public $cref/sample_random/base_class/sample_random/$$
$srccode%cpp% */
	std::string try_sample_random(
		d_vector&             sample               ,
		const std::string&    random_ipopt_options ,
		const d_vector&       fixed_vec            ,
		const d_vector&       random_lower         ,
		const d_vector&       random_upper         ,
		const d_vector&       random_in
	);
/* %$$
$comment ------------------------------------------------------------------- $$
$head Other Member Functions$$

$subhead fix_con_eval$$
See $cref fix_con_eval$$
$srccode%cpp% */
	d_vector fix_con_eval(const d_vector& fixed_vec);
	friend bool ::fix_con_eval_xam(void);
/* %$$

$subhead fix_con_jac$$
See $cref fix_con_jac$$
$srccode%cpp% */
	void fix_con_jac(
		const d_vector&        fixed_vec   ,
		CppAD::vector<size_t>& row_out     ,
		CppAD::vector<size_t>& col_out     ,
		d_vector&              val_out
	);
	friend bool ::fix_con_jac_xam(void);
/* %$$

$subhead fix_con_hes$$
See $cref fix_con_hes$$
$srccode%cpp% */
	void fix_con_hes(
		const d_vector&        fixed_vec   ,
		const d_vector&        weight      ,
		CppAD::vector<size_t>& row_out     ,
		CppAD::vector<size_t>& col_out     ,
		d_vector&              val_out
	);
	friend bool ::fix_con_hes_xam(void);
/* %$$

$subhead fix_like_eval$$
See $cref fix_like_eval$$
$srccode%cpp% */
	d_vector fix_like_eval(const d_vector& fixed_vec);
	friend bool ::fix_like_eval_xam(void);
/* %$$

$subhead fix_like_jac$$
See $cref fix_like_jac$$
$srccode%cpp% */
	void fix_like_jac(
		const d_vector&        fixed_vec   ,
		CppAD::vector<size_t>& row_out     ,
		CppAD::vector<size_t>& col_out     ,
		d_vector&              val_out
	);
	friend bool ::fix_like_jac_xam(void);
/* %$$

$subhead fix_like_hes$$
See $cref fix_like_hes$$
$srccode%cpp% */
	void fix_like_hes(
		const d_vector&        fixed_vec   ,
		const d_vector&        weight      ,
		CppAD::vector<size_t>& row_out     ,
		CppAD::vector<size_t>& col_out     ,
		d_vector&              val_out
	);
	friend bool ::fix_like_hes_xam(void);
/* %$$

$subhead logdet_jac$$
See $cref logdet_jac$$
$srccode%cpp% */
	void logdet_jac(
		const d_vector& fixed_vec  ,
		const d_vector& random_vec ,
		d_vector&       logdet_fix ,
		d_vector&       logdet_ran
	);
	friend bool ::logdet_jac_xam(void);
/* %$$

$subhead ran_con_eval$$
See $cref ran_con_eval$$
$srccode%cpp% */
	void ran_con_eval(
		const d_vector& random_vec ,
		d_vector&       Au
	);
	friend bool ::ran_con_eval_xam(void);
/* %$$

$subhead ran_con_jac$$
See $cref ran_con_jac$$
$srccode%cpp% */
	void ran_con_jac(
		const d_vector&                fixed_vec  ,
		const d_vector&                random_vec ,
		d_sparse_rcv&                  jac_info
	);
	friend bool ::ran_con_jac_xam(void);
	friend bool ::sample_fixed_1(void);
/* %$$

$subhead ran_obj_eval$$
See $cref ran_obj_eval$$
$srccode%cpp% */
	double ran_obj_eval(
		const d_vector& fixed_vec  ,
		const d_vector& random_vec
	);
	friend bool ::ran_obj_eval_xam(void);
	friend bool ::laplace_obj_fun(void);
	friend bool ::laplace_obj_tst(void);
/* %$$

$subhead ran_obj_jac$$
See $cref ran_obj_jac$$
$srccode%cpp% */
	void ran_obj_jac(
		const d_vector& fixed_vec  ,
		const d_vector& random_vec ,
		d_vector&       r_fixed
	);
	friend bool ::ran_obj_jac_xam(void);
	friend bool ::der_var_hes(void);
	friend bool ::delta_ran_obj(void);
/* %$$

$subhead laplace_obj_hes$$
See $cref laplace_obj_hes$$
$srccode%cpp% */
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
/* %$$

$subhead update_factor$$
See $cref update_factor$$
$srccode%cpp% */
	void update_factor(
		const d_vector&         fixed_vec   ,
		const d_vector&         random_vec
	);
	friend bool ::update_factor_xam(void);
/* %$$
$end
-------------------------------------------------------------------------------
*/
# endif
