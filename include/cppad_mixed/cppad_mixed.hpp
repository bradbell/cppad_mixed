// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_CPPAD_MIXED_HPP
# define CPPAD_MIXED_CPPAD_MIXED_HPP

# include <map>
# include <cppad/cppad.hpp>
# include <cppad_mixed/newton_step.hpp>
# include <cppad_mixed/sparse_hes_info.hpp>
# include <cppad_mixed/sparse_jac_info.hpp>

// private examples
extern bool constraint_eval_xam(void);
extern bool constraint_jac_xam(void);
extern bool constraint_hes_xam(void);
extern bool ran_like_grad_xam(void);
extern bool ranobj_hes_xam(void);
extern bool fix_like_eval_xam(void);
extern bool fix_like_jac_xam(void);
extern bool fix_like_hes_xam(void);
extern bool hes_ran_fun_xam(void);
extern bool ranobj_eval_xam(void);
extern bool logdet_grad_xam(void);
extern bool ranobj_grad_xam(void);
extern bool hes_cross_xam(void);

//  tests
extern bool der_var_hes(void);
extern bool delta_ranobj(void);

namespace cppad_mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

class cppad_mixed {
	friend class optimize_random_eval;
	friend class ipopt_fixed;
public:
/*
$begin cppad_mixed_public$$
$spell
	init
	ranobj
	cppad
	obj
	const
	vec
	typedef
	CppAD
	ctor
	std
	bool
$$

$section cppad_mixed Public Declarations$$

$head AD Types$$
$index a_double$$
$codep */
	typedef CppAD::AD<double>          a1_double;
	typedef CppAD::AD<a1_double>       a2_double;
/* $$

$head Vector Types$$
$mindex d_vector a1d_vector a2d_vector$$
$codep */
	typedef CppAD::vector<double>      d_vector;
	typedef CppAD::vector<a1_double>   a1d_vector;
	typedef CppAD::vector<a2_double>   a2d_vector;
/* $$
$head User Defined$$
The following are $code cppad_mixed$$ pure virtual functions and hence must
be defined by the user's derived class:

$subhead ran_like$$
$codep */
	virtual a2d_vector ran_like(
		const a2d_vector& fixed_vec  ,
		const a2d_vector& random_vec
	) = 0;
	virtual a1d_vector ran_like(
		const a1d_vector& fixed_vec  ,
		const a1d_vector& random_vec
	) = 0;
/* $$
See $cref/ran_like/cppad_mixed_ran_like/$$.

$subhead fix_like$$
$codep */
	virtual a1d_vector fix_like(
		const a1d_vector& fixed_vec
	) = 0 ;
/* $$
See $cref/fix_like/cppad_mixed_fix_like/$$.

$subhead constraint$$
$codep */
	virtual a1d_vector constraint(
		const a1d_vector& fixed_vec
	) = 0 ;
/* $$
See $cref/constraint/cppad_mixed_constraint/$$.

$subhead fatal_error$$
$codep */
	virtual void fatal_error(const std::string& error_message) = 0;
/* $$
This routine displays error message and then exits the program
(hence it does not return).

$subhead warning$$
$codep */
	virtual void warning(const std::string& warning_message) = 0;
/* $$
This routine displays a warning message and then returns.

$head constructor$$
Construct an $code cppad_mixed$$ derived class object; see
$cref/derived_ctor/cppad_mixed_derived_ctor/$$.
$codep */
	cppad_mixed(size_t n_fixed, size_t n_random, bool quasi_fixed)
/* $$
$comment */
	:
	n_fixed_(n_fixed)               ,
	n_random_(n_random)             ,
	quasi_fixed_(quasi_fixed)       ,
	initialize_done_(false)         ,
	init_fix_like_done_(false)    ,
	init_constraint_done_(false)  ,
	init_ran_like_done_(false)    ,
	init_hes_ran_done_(false)     ,
	init_hes_cross_done_(false)   ,
	record_newton_atom_done_(false) ,
	init_ranobj_done_(false)     ,
	init_hes_ranobj_done_(false)
	{ }
/* $$
$head initialize$$
Directly after construction, use this function to initialize
the derived class object; see $cref/initialize/cppad_mixed_initialize/$$.
$codep */
	std::map<std::string, size_t> initialize(
		const d_vector& fixed_vec, const d_vector& random_vec
	);
/* $$
$head optimize_random$$
Given the fixed effects, optimize with respect to the random effects;
see  $cref/optimize_random/cppad_mixed_optimize_random/$$.
$codep */
	d_vector optimize_random(
		const std::string& options      ,
		const d_vector&    fixed_vec    ,
		const d_vector&    random_lower ,
		const d_vector&    random_upper ,
		const d_vector&    random_in
	);
/* $$
$head optimize_fixed$$
Optimize the total objective with respect to the fixed effects;
see  $cref/optimize_fixed/cppad_mixed_optimize_fixed/$$.
$codep */
	d_vector optimize_fixed(
		const std::string& fixed_options    ,
		const std::string& random_options   ,
		const d_vector&    fixed_lower      ,
		const d_vector&    fixed_upper      ,
		const d_vector&    constraint_lower ,
		const d_vector&    constraint_upper ,
		const d_vector&    fixed_in         ,
		const d_vector&    random_lower     ,
		const d_vector&    random_upper     ,
		const d_vector&    random_in
	);
/* $$
$childtable%
	src/derived_ctor.omh%
	src/initialize.cpp%
	src/ran_like.omh%
	src/fix_like.omh%
	src/constraint.omh%
	src/optimize_random.cpp%
	src/optimize_fixed.cpp
%$$
$end
*/
private:
/*
------------------------------------------------------------------------------
$begin cppad_mixed_private$$
$spell
	init
	ranobj
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
$$

$section cppad_mixed Private Declarations$$

$childtable%include/cppad_mixed/mixed_pack.hpp
	%include/cppad_mixed/mixed_unpack.hpp
	%src/constraint_eval.cpp
	%src/constraint_jac.cpp
	%src/constraint_hes.cpp
	%src/init_ran_like.cpp
	%src/init_hes_ran.cpp
	%src/init_hes_cross.cpp
	%src/newton_step.cpp
	%src/init_ranobj.cpp
	%src/init_hes_ranobj.cpp
	%src/init_fix_like.cpp
	%src/init_constraint.cpp
	%src/ranobj_eval.cpp
	%src/logdet_grad.cpp
	%src/ranobj_grad.cpp
	%src/ran_like_grad.cpp
	%src/ranobj_hes.cpp
	%src/fix_like_eval.cpp
	%src/fix_like_jac.cpp
	%src/fix_like_hes.cpp
	%src/chol_hes_ran.cpp
	%include/cppad_mixed/sparse_hes_info.hpp
	%include/cppad_mixed/sparse_jac_info.hpp
%$$

$head n_fixed_$$
The number of fixed effects is given by
$codep */
	const size_t n_fixed_;
/* $$
$head n_random_$$
The number of random effects is given by
$codep */
	const size_t n_random_;
/* $$
$head quasi_fixed_$$
Are we using a quasi-Newton method (or full Newton method)
when $cref/optimizing fixed effects/cppad_mixed_optimize_fixed/$$.
$codep */
	const bool quasi_fixed_;
/* $$
$head initialize_done_$$
The following flag is false after construction and true after
the corresponding member function is called:
$codep */
	bool                initialize_done_;
	bool                init_fix_like_done_;
	bool                init_constraint_done_;
	// only called when n_random_ > 0
	bool                init_ran_like_done_;
	bool                init_hes_ran_done_;
	bool                init_hes_cross_done_;
	// only called when n_random_ > 0 and quasi_fixed_ is false
	bool                record_newton_atom_done_;
	bool                init_ranobj_done_;
	bool                init_hes_ranobj_done_;
/* $$

$head ran_like$$
If $icode%n_random_% > 0%$$ and $code init_ran_like_done_$$,
$cref/ran_like_fun_/init_ran_like/ran_like_fun_/$$ and
$cref/ran_like_a1fun_/init_ran_like/ran_like_a1fun_/$$ are
recordings of the user's $cref/ran_like/cppad_mixed_ran_like/$$.
function.
$codep */
	CppAD::ADFun<double>      ran_like_fun_;
	CppAD::ADFun<a1_double>   ran_like_a1fun_;
/* $$
The following objects hold information for computing derivatives
with these ADFun objects:

$subhead hes_ran_$$
If $icode%n_random_% > 0%$$ and $code init_hes_ran_done_$$,
$cref/hes_ran_/init_hes_ran/hes_ran_/$$ contains
information for the Hessian of the
$cref/random likelihood
	/cppad_mixed_theory
	/Random Likelihood, f(theta, u)
/$$
with respect to the random effects; i.e.
$latex f_{uu}^{(2)} ( \theta , u )$$.
$codep */
	sparse_hes_info             hes_ran_;
	// recording of sparse Hessian calculation
	CppAD::ADFun<double>        hes_ran_fun_;
	//
	friend bool ::hes_ran_fun_xam(void);
/* $$

$subhead hes_cross_$$
If $icode%n_random_% > 0%$$ and $code init_hes_cross_done_$$,
$cref/hes_cross_/init_hes_cross/hes_cross_/$$ contains
information for the cross partials of the Hessian of the
$cref/random likelihood
	/cppad_mixed_theory
	/Random Likelihood, f(theta, u)
/$$
; i.e.  $latex f_{u \theta}^{(2)} ( \theta , u )$$.
$codep */
	sparse_hes_info             hes_cross_;
	//
	friend bool ::hes_cross_xam(void);
/* $$


$head newton_atom_$$
If $icode%n_random_% > 0%$$, quasi_fixed_ is false, and
$code record_newton_atom_done_$$,
this is a CppAD atomic function that computes one Newton Step in the
solution of the equation $latex f_u ( \theta, u) = 0$$ as well
as the log of the determinant of $latex f_{uu} ( \theta , u )$$;
see $cref/initialize newton_step/newton_step/initialize/$$.
$codep */
	// computation of the Hessian as an atomic operation
	newton_step                 newton_atom_;
/* $$

$head ranobj_fun_$$
If $icode%n_random_% > 0%$$, quasi_fixed_ is false, and
$code init_ranobj_done_$$,
this is a recording of the second approximation for the
random part of the Laplace approximation, $latex H( \beta , \theta , u)$$;
see $cref/ranobj_fun_/init_ranobj/ranobj_fun_/$$.
$codep */
	CppAD::ADFun<double>        ranobj_fun_;   // for computing H_beta_beta
/* $$
The following objects hold information for computing derivatives
with this ADFun object:

$subhead hes_ranobj_$$
If $icode%n_random_% > 0%$$, quasi_fixed_ is false, and
$code init_hes_ranobj_done_$$,
$cref/hes_ranobj_/init_hes_ranobj/hes_ranobj_/$$ contains
information for the Hessian of the
$cref/random objective
	/cppad_mixed_theory
	/Objective
	/Random Objective, r(theta)
/$$
$codep */
	sparse_hes_info             hes_ranobj_;
/* $$

$head fix_like_$$
$cref/fix_like_fun_/init_fix_like/fix_like_fun_/$$
is a recording of the fixed part of the likelihood function; see,
$cref/fix_like/cppad_mixed_fix_like/$$.
$codep */
	CppAD::ADFun<double>        fix_like_fun_;     // g(theta)
/* $$
The following objects hold information for computing derivatives
with this ADFun object:

$subhead fix_like_jac_$$
$cref/fix_like_jac_/init_fix_like/fix_like_jac_/$$
contains information for the Jacobian of the
$cref/fixed likelihood/cppad_mixed_theory/Fixed Likelihood, g(theta)/$$.
$codep */
	sparse_jac_info             fix_like_jac_;
/* $$

$subhead fix_like_hes_$$
If $icode quasi_fixed$$ is false,
$cref/fix_like_hes_/init_fix_like/fix_like_hes_/$$
contains information for the Hessian of the
$cref/fixed likelihood/cppad_mixed_theory/Fixed Likelihood, g(theta)/$$.
$codep */
	sparse_hes_info             fix_like_hes_;
/* $$

$head constraint_fun_$$
$cref/constraint_fun_/init_constraint/constraint_fun_/$$
is a recording of the fixed part of the likelihood function; see,
$cref/constraint/cppad_mixed_constraint/$$.
$codep */
	CppAD::ADFun<double>        constraint_fun_;     // c(theta)
/* $$
The following objects hold information for computing derivatives
with this ADFun object:

$subhead constraint_jac_$$
$cref/constraint_jac_/init_constraint/constraint_jac_/$$
contains information for the Jacobian of the
constraint function $latex c ( \theta )$$.
$codep */
	sparse_jac_info             constraint_jac_;
/* $$

$subhead constraint_hes_$$
If $icode quasi_fixed$$ is false,
$cref/constraint_hes_/init_constraint/constraint_hes_/$$
contains information for the Hessian of the
$cref/constraints/cppad_mixed_constraint/$$ function $latex c( \theta )$$.
The corresponding ADFun object is
$cref/constraint_fun_/init_constraint/constraint_fun_/$$.
$codep */
	sparse_hes_info             constraint_hes_;
/* $$
------------------------------------------------------------------------------
$head pack$$
See $cref cppad_mixed_pack$$.
$codep */
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
/* $$
$head unpack$$
See $cref cppad_mixed_unpack$$.
$codep */
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
/* $$
------------------------------------------------------------------------------
$head init_ran_like$$
See $cref init_ran_like$$.
$codep */
	void init_ran_like(
		const d_vector& fixed_vec ,
		const d_vector& random_vec
	);
/* $$
$head init_hes_ran$$
See $cref init_hes_ran$$.
$codep */
	void init_hes_ran(
		const d_vector& fixed_vec ,
		const d_vector& random_vec
	);
/* $$
$head init_hes_cross$$
See $cref init_hes_cross$$.
$codep */
	void init_hes_cross(
		const d_vector& fixed_vec ,
		const d_vector& random_vec
	);
/* $$
$head init_ranobj$$
See $cref init_ranobj$$.
$codep */
	void init_ranobj(
		const d_vector& fixed_vec ,
		const d_vector& random_vec
	);
/* $$
$head init_hes_ranobj$$
See $cref init_hes_ranobj$$.
$codep */
	void init_hes_ranobj(
		const d_vector& fixed_vec ,
		const d_vector& random_vec
	);
/* $$
$head init_fix_like$$
See $cref init_fix_like$$.
$codep */
	void init_fix_like(const d_vector& fixed_vec);
/* $$
$head init_constraint$$
See $cref init_constraint$$.
$codep */
	void init_constraint(const d_vector& fixed_vec);
/* $$
------------------------------------------------------------------------------
$head ranobj_eval$$
See $cref ranobj_eval$$
$codep */
	double ranobj_eval(
		const d_vector& fixed_vec  ,
		const d_vector& random_vec
	);
	friend bool ::ranobj_eval_xam(void);
/* $$
------------------------------------------------------------------------------
$head logdet_grad$$
See $cref logdet_grad$$
$codep */
	void logdet_grad(
		const d_vector& fixed_vec  ,
		const d_vector& random_vec ,
		d_vector&       logdet_fix ,
		d_vector&       logdet_ran
	);
	friend bool ::logdet_grad_xam(void);
/* $$
------------------------------------------------------------------------------
$head ranobj_grad$$
See $cref ranobj_grad$$
$codep */
	void ranobj_grad(
		const d_vector& fixed_vec  ,
		const d_vector& random_vec ,
		d_vector&       r_fixed
	);
	friend bool ::ranobj_grad_xam(void);
	friend bool ::der_var_hes(void);
	friend bool ::delta_ranobj(void);
/* $$
------------------------------------------------------------------------------
$head ran_like_grad$$
See $cref ran_like_grad$$
$codep */
	// ran_like_grad
	a1d_vector ran_like_grad(
		const a1d_vector&       fixed_vec   ,
		const a1d_vector&       random_vec
	);
	friend bool ::ran_like_grad_xam(void);
/* $$
------------------------------------------------------------------------------
$head ranobj_hes$$
See $cref ranobj_hes$$
$codep */
	// ranobj_hes
	void ranobj_hes(
		const d_vector&         fixed_vec   ,
		const d_vector&         random_vec  ,
		CppAD::vector<size_t>&  row_out     ,
		CppAD::vector<size_t>&  col_out     ,
		d_vector&               val_out
	);
	friend bool ::ranobj_hes_xam(void);
/* $$
-------------------------------------------------------------------------------
$head constraint_eval$$
See $cref constraint_eval$$
$codep */
	// constraint_eval
	d_vector constraint_eval(const d_vector& fixed_vec);
	friend bool ::constraint_eval_xam(void);
/* $$
$end
$head constraint_jac$$
See $cref constraint_jac$$
$codep */
	// constraint_jac
	void constraint_jac(
		const d_vector&        fixed_vec   ,
		CppAD::vector<size_t>& row_out     ,
		CppAD::vector<size_t>& col_out     ,
		d_vector&              val_out
	);
	friend bool ::constraint_jac_xam(void);
/* $$
$end
$head constraint_hes$$
See $cref constraint_hes$$
$codep */
	// constraint_hes
	void constraint_hes(
		const d_vector&        fixed_vec   ,
		const d_vector&        weight      ,
		CppAD::vector<size_t>& row_out     ,
		CppAD::vector<size_t>& col_out     ,
		d_vector&              val_out
	);
	friend bool ::constraint_hes_xam(void);
/* $$
-------------------------------------------------------------------------------
$head fix_like_eval$$
See $cref fix_like_eval$$
$codep */
	// fix_like_eval
	d_vector fix_like_eval(const d_vector& fixed_vec);
	friend bool ::fix_like_eval_xam(void);
/* $$
$end
$head fix_like_jac$$
See $cref fix_like_jac$$
$codep */
	// fix_like_jac
	void fix_like_jac(
		const d_vector&        fixed_vec   ,
		CppAD::vector<size_t>& row_out     ,
		CppAD::vector<size_t>& col_out     ,
		d_vector&              val_out
	);
	friend bool ::fix_like_jac_xam(void);
/* $$
$end
$head fix_like_hes$$
See $cref fix_like_hes$$
$codep */
	// fix_like_hes
	void fix_like_hes(
		const d_vector&        fixed_vec   ,
		const d_vector&        weight      ,
		CppAD::vector<size_t>& row_out     ,
		CppAD::vector<size_t>& col_out     ,
		d_vector&              val_out
	);
	friend bool ::fix_like_hes_xam(void);
/* $$
$end
-------------------------------------------------------------------------------
*/
};

} // END_CPPAD_MIXED_NAMESPACE

# include <cppad_mixed/mixed_pack.hpp>
# include <cppad_mixed/mixed_unpack.hpp>

# endif
