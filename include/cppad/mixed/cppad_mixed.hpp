// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_CPPAD_MIXED_HPP
# define CPPAD_MIXED_CPPAD_MIXED_HPP

# include <map>
# include <cppad/cppad.hpp>
# include <cppad/mixed/sparse_hes_rcv.hpp>
# include <cppad/mixed/sparse_jac_rcv.hpp>
# include <cppad/mixed/ldlt_eigen.hpp>
# include <cppad/mixed/sparse_mat_info.hpp>
# include <cppad/mixed/configure.hpp>
# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cppad/mixed/ldlt_eigen.hpp>
# include <cppad/mixed/fixed_solution.hpp>
# include <cppad/mixed/typedef.hpp>

/*
Examples and tests that use cppad_mixed private information.
This list can be generated using the following command:
grep -P '^\t*friend *bool *::' cppad_mixed.hpp | \
	sed -e 's|^\t*friend *bool *::|extern bool |' | sort
*/
extern bool delta_ran_obj(void);
extern bool der_var_hes(void);
extern bool fix_con_eval_xam(void);
extern bool fix_con_hes_xam(void);
extern bool fix_con_jac_xam(void);
extern bool fix_like_eval_xam(void);
extern bool fix_like_hes_xam(void);
extern bool fix_like_jac_xam(void);
extern bool hes_cross_xam(void);
extern bool laplace_obj_tst(void);
extern bool laplace_obj_hes(void);
extern bool laplace_obj_hes_xam(void);
extern bool logdet_jac_xam(void);
extern bool ran_con_eval_xam(void);
extern bool ran_con_jac_xam(void);
extern bool ran_hes_fun_xam(void);
extern bool ran_jac_fun_xam(void);
extern bool ran_like_jac_xam(void);
extern bool ran_obj_eval_xam(void);
extern bool ran_obj_jac_xam(void);
extern bool ran_obj_tst(void);
extern bool sample_fixed_1(void);
extern bool update_factor_xam(void);


namespace CppAD { namespace mixed {
	class optimize_random_ipopt;
	class ipopt_fixed;
	class ipopt_random;
} }

class cppad_mixed {
	friend class CppAD::mixed::optimize_random_ipopt;
	friend class CppAD::mixed::ipopt_fixed;
	friend class CppAD::mixed::ipopt_random;
	friend bool ::ran_obj_tst(void);
public:
/*
$begin public$$
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
	typedef CppAD::mixed::a2_double      a2_double;
	typedef CppAD::mixed::a3_double      a3_double;
	// vector types
	typedef CppAD::mixed::s_vector       s_vector;
	typedef CppAD::mixed::d_vector       d_vector;
	typedef CppAD::mixed::a1_vector      a1_vector;
	typedef CppAD::mixed::a2_vector      a2_vector;
	typedef CppAD::mixed::a3_vector      a3_vector;
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
	virtual a3_vector ran_likelihood(
		const a3_vector& fixed_vec  ,
		const a3_vector& random_vec )
	{	return a3_vector(0); }
/* %$$
See $cref ran_likelihood$$.

$subhead ran_likelihood_jac$$
This function is only used to increase speed and reduce memory use.
$srccode%cpp% */
	virtual a2_vector ran_likelihood_jac(
		const a2_vector&           fixed_vec  ,
		const a2_vector&           random_vec )
	{	return a2_vector(0); }
/* %$$
See $cref ran_likelihood_jac$$.

$subhead ran_likelihood_hes$$
This function is only used to increase speed and reduce memory use.
$srccode%cpp% */
	virtual a1_vector ran_likelihood_hes(
		const a1_vector&             fixed_vec  ,
		const a1_vector&             random_vec ,
		const CppAD::vector<size_t>& row        ,
		const CppAD::vector<size_t>& col        )
	{	return a1_vector(0); }
/* %$$
See $cref ran_likelihood_hes$$.

$subhead fix_likelihood$$
This function should be used if there is a prior on the fixed effects,
or there is data that does not depend on the random effects.
$srccode%cpp% */
	virtual a1_vector fix_likelihood(
		const a1_vector& fixed_vec )
	{	return a1_vector(0); }
/* %$$
See $cref fix_likelihood$$.

$subhead fix_constraint$$
This function is used to define constraints
that only depend on the fixed effects.
$srccode%cpp% */
	virtual a1_vector fix_constraint(
		const a1_vector& fixed_vec )
	{	return a1_vector(0); }
/* %$$
See $cref fix_constraint$$.

$subhead fatal_error$$
This routine displays an error message and then exits the program.
Its default definition below can be replaced by a user definition.
Note that if $code NDEBUG$$ is not defined, this generates an assert,
otherwise it exits.
$srccode%cpp% */
	virtual void fatal_error(const std::string& error_message)
	{	std::cerr << "cppad_mixed error: " << error_message << std::endl;
		assert(false);
		exit(1);
	}

/* %$$

$subhead warning$$
This routine displays a warning message and then returns.
Its default definition below can be replaced by a user definition.
$srccode%cpp% */
	virtual void warning(const std::string& warning_message)
	{	std::cerr << "cppad_mixed warning: " << warning_message << std::endl;
	}
/* %$$

$head constructor$$
$cref derived_ctor$$, $title derived_ctor$$.
$srccode%cpp% */
	cppad_mixed(
		size_t               n_fixed       ,
		size_t               n_random      ,
		bool                 quasi_fixed   ,
		bool                 bool_sparsity ,
		const d_sparse_rcv&  A_rcv
	);
/* %$$

$head destructor$$
$srccode%cpp% */
	~cppad_mixed(void);
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
		const d_vector&    random_in
	);
/* %$$
$head information_mat$$
$cref information_mat$$, $title information_mat$$.
$srccode%cpp% */
	d_sparse_rcv information_mat(
		const CppAD::mixed::fixed_solution&  solution             ,
		const d_vector&                      random_opt
	);
/* %$$
$head sample_fixed$$
$cref sample_fixed$$, $title sample_fixed$$.
$srccode%cpp% */
	void sample_fixed(
		d_vector&                            sample               ,
		const d_sparse_rcv&                  information_rcv      ,
		const CppAD::mixed::fixed_solution&  solution             ,
		const d_vector&                      fixed_lower          ,
		const d_vector&                      fixed_upper          ,
		const d_vector&                      random_opt
	);
/* %$$
$head sample_random$$
$cref sample_random$$, $title sample_random$$.
$srccode%cpp% */
	void sample_random(
		d_vector&            sample               ,
		const std::string&   random_ipopt_options ,
		const d_vector&      fixed_vec            ,
		const d_vector&      random_lower         ,
		const d_vector&      random_upper         ,
		const d_vector&      random_in
	);
/* %$$
$childtable%src/derived_ctor.cpp
	%src/ran_likelihood.omh
	%src/fix_likelihood.omh
	%src/fix_constraint.omh
	%src/initialize.cpp
	%src/optimize_random.cpp
	%src/optimize_fixed.cpp
	%src/eigen/information_mat.cpp
	%src/eigen/sample_fixed.cpp
	%src/sample_random.cpp
	%src/ran_likelihood_jac.omh
	%src/ran_likelihood_hes.omh
%$$


$end
*/
private:
/*
------------------------------------------------------------------------------
$begin private$$
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

	%src/init_ran_jac.cpp
	%src/init_ran_hes.cpp
	%src/init_laplace_obj.cpp
	%src/init_ldlt_ran_hes.cpp
	%src/init_fix_con.cpp
	%src/init_fix_like.cpp
	%src/init_hes_cross.cpp
	%src/init_laplace_obj_hes.cpp
	%src/init_ran_like.cpp

	%src/fix_con_eval.cpp
	%src/fix_con_hes.cpp
	%src/fix_con_jac.cpp
	%src/fix_like_eval.cpp
	%src/fix_like_hes.cpp
	%src/fix_like_jac.cpp
	%src/logdet_jac.cpp
	%src/ran_like_jac.cpp
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
	bool                init_newton_checkpoint_done_;
	bool                init_laplace_obj_done_;
	bool                init_laplace_obj_hes_done_;
	// called in all cases
	bool                init_fix_like_done_;
	bool                init_fix_con_done_;
	// true when all initialization (for this case) is done
	bool                initialize_done_;

/* %$$
$head cppad_error_handler$$
Used to map CppAD error messages to
$cref/fatal_error/public/User Defined Functions/fatal_error/$$.
$srccode%cpp% */
	CppAD::ErrorHandler cppad_error_handler_;
/* %$$

$head ran_like_fun_$$
If $icode%n_random_% > 0%$$ and $code init_ran_like_done_$$,
$cref/ran_like_fun_/init_ran_like/ran_like_fun_/$$,
$cref/ran_like_a1fun_/init_ran_like/ran_like_a1fun_/$$, and
$cref/ran_like_a2fun_/init_ran_like/ran_like_a2fun_/$$,
are recordings of the user's $cref ran_likelihood$$.
function.
$srccode%cpp% */
	CppAD::ADFun<double>      ran_like_fun_;
	CppAD::ADFun<a1_double>   ran_like_a1fun_;
	CppAD::ADFun<a2_double>   ran_like_a2fun_;
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
	CppAD::ADFun<a1_double>  ran_jac_a1fun_;
	sparse_rc                ran_jac2hes_rc_;
	//
	friend bool ::ran_jac_fun_xam(void);
/* %$$
The row indices in $icode ran_jac2hes_rc_$$
are for just the random effects $latex u$$,
the column indices are for both fixed and random effects
$latex ( \theta , u )$$.

$head ran_hes_fun_$$
If $icode%n_random_% > 0%$$ and $code init_ran_hes_done_$$,
$cref/ran_hes_rcv_/init_ran_hes/ran_hes_rcv_/$$
contains information for the Hessian of the
$cref/random likelihood
	/theory
	/Random Likelihood, f(theta, u)
/$$
with respect to the random effects; i.e.
$latex f_{u,u} ( \theta , u )$$.
$srccode%cpp% */
	// sparse Hessian
	d_sparse_rcv                  ran_hes_rcv_;
	//
	// recording of sparse Hessian calculation
	CppAD::ADFun<double>        ran_hes_fun_;
	//
	friend bool ::ran_hes_fun_xam(void);
/* %$$
The row and column indices in $icode ran_hes_rcv_$$
are for both the fixed and random effects $latex ( \theta , u )$$
and hence are all greater than or equal  $icode n_fixed$$.

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
	CPPAD_MIXED_LDLT_CLASS ldlt_ran_hes_;
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
$code init_laplace_obj_done_$$,
this is a recording of the second order approximation for the
random part of the Laplace approximation, $latex H( \beta , \theta , u)$$;
see $cref/laplace_obj_fun_/init_laplace_obj/laplace_obj_fun_/$$.
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

$subhead check_user_ran_hes$$
See $cref check_user_ran_hes$$.
$srccode%cpp% */
	void check_user_ran_hes(
		const d_vector& fixed_vec ,
		const d_vector& random_vec
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

$subhead init_laplace_obj$$
See $cref init_laplace_obj$$.
$srccode%cpp% */
	void init_laplace_obj(
		const d_vector& fixed_vec ,
		const d_vector& random_vec
	);
/* %$$
$comment ------------------------------------------------------------------- $$
$head Try Member Functions$$
For each $codei%try_%name%$$ case below,
the $cref public$$ function $icode name$$ calls the private function
$codei%try_%name%$$ from a $code try$$ block
with a corresponding $code catch$$ that maps a
$code cppad_mixed$$ $cref exception$$ to a
$cref/fatal_error/public/User Defined Functions/fatal_error/$$ call.

$subhead try_initialize$$
Called by public $cref/initialize/public/initialize/$$
$srccode%cpp% */
	std::map<std::string, size_t> try_initialize(
		const d_vector&  fixed_vec  ,
		const d_vector&  random_vec
	);
/* %$$
$subhead try_optimize_random$$
Called by public $cref/optimize_random/public/optimize_random/$$
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
Called by public $cref/optimize_fixed/public/optimize_fixed/$$
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
$subhead try_information_mat$$
Called by public $cref/information_mat/public/information_mat/$$
$srccode%cpp% */
	d_sparse_rcv try_information_mat(
		const CppAD::mixed::fixed_solution&  solution      ,
		const d_vector&                      random_opt
	);
/* %$$
$subhead try_sample_fixed$$
Called by public $cref/sample_fixed/public/sample_fixed/$$
$srccode%cpp% */
	void try_sample_fixed(
		d_vector&                            sample               ,
		const d_sparse_rcv&                  information_rcv      ,
		const CppAD::mixed::fixed_solution&  solution             ,
		const d_vector&                      fixed_lower          ,
		const d_vector&                      fixed_upper          ,
		const d_vector&                      random_opt
	);
/* %$$
$subhead sample_random$$
Called by public $cref/sample_random/public/sample_random/$$
$srccode%cpp% */
	void try_sample_random(
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

$subhead ran_like_jac$$
See $cref ran_like_jac$$
$srccode%cpp% */
	a1_vector ran_like_jac(
		const a1_vector&        fixed_vec   ,
		const a1_vector&        random_vec
	);
	friend bool ::ran_like_jac_xam(void);
/* %$$

$subhead check_user_ran_jac$$
See $cref check_user_ran_jac$$
$srccode%cpp% */
	void check_user_ran_jac(
		const d_vector&       fixed_vec     ,
		const d_vector&       random_vec
	);
/* %$$

$subhead ran_obj_eval$$
See $cref ran_obj_eval$$
$srccode%cpp% */
	double ran_obj_eval(
		const d_vector& fixed_vec  ,
		const d_vector& random_vec
	);
	friend bool ::ran_obj_eval_xam(void);
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
};


# include <cppad/mixed/pack.hpp>
# include <cppad/mixed/unpack.hpp>

# endif
