// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_CPPAD_MIXED_HPP
# define CPPAD_MIXED_CPPAD_MIXED_HPP

# include <map>
# include <cppad/cppad.hpp>
# include <cppad/mixed/newton_step.hpp>
# include <cppad/mixed/sparse_hes_info.hpp>
# include <cppad/mixed/sparse_jac_info.hpp>

// private examples
extern bool fix_con_eval_xam(void);
extern bool fix_con_hes_xam(void);
extern bool fix_con_jac_xam(void);
extern bool fix_like_eval_xam(void);
extern bool fix_like_hes_xam(void);
extern bool fix_like_jac_xam(void);
extern bool hes_cross_xam(void);
extern bool hes_ran_fun_xam(void);
extern bool logdet_jac_xam(void);
extern bool ran_like_jac_xam(void);
extern bool ran_obj_eval_xam(void);
extern bool ran_obj_jac_xam(void);
extern bool ran_obj_hes_xam(void);

//  tests
extern bool der_var_hes(void);
extern bool delta_ran_obj(void);


namespace CppAD { namespace mixed {
	class optimize_random_eval;
	class ipopt_fixed;
} }

class cppad_mixed {
	friend class CppAD::mixed::optimize_random_eval;
	friend class CppAD::mixed::ipopt_fixed;
public:
/*
$begin public$$
$spell
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
$$

$section cppad_mixed: Public Declarations$$
These $code cppad_mixed$$ class declarations are $code public$$.
They are part of the user API and can be used by a derived class object
$cref/mixed_object/derived_ctor/mixed_object/$$.

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

$subhead ran_likelihood$$
$codep */
	virtual a2d_vector ran_likelihood(
		const a2d_vector& fixed_vec  ,
		const a2d_vector& random_vec )
	{	return a2d_vector(0); }
	virtual a1d_vector ran_likelihood(
		const a1d_vector& fixed_vec  ,
		const a1d_vector& random_vec )
	{	return a1d_vector(0); }
/* $$
See $cref ran_likelihood$$.

$subhead fix_likelihood$$
$codep */
	virtual a1d_vector fix_likelihood(
		const a1d_vector& fixed_vec )
	{	return a1d_vector(0); }
/* $$
See $cref fix_likelihood$$.

$subhead fix_constraint$$
$codep */
	virtual a1d_vector fix_constraint(
		const a1d_vector& fixed_vec )
	{	return a1d_vector(0); }
/* $$
See $cref fix_constraint$$.

$subhead fatal_error$$
This routine displays an error message and then exits the program.
Its default definition below can be replaced by a user definition.
$codep */
	virtual void fatal_error(const std::string& error_message)
	{	std::cerr << "cppad_mixed error: " << error_message << std::endl;
		assert(false);
	}

/* $$

$subhead warning$$
This routine displays a warning message and then returns.
Its default definition below can be replaced by a user definition.
$codep */
	virtual void warning(const std::string& warning_message)
	{	std::cerr << "cppad_mixed warning: " << warning_message << std::endl;
		assert(false);
	}
/* $$

$head constructor$$
Construct an $code cppad_mixed$$ derived class object; see
$cref/derived_ctor/derived_ctor/$$.
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
	init_fix_con_done_(false)  ,
	init_ran_like_done_(false)    ,
	init_hes_ran_done_(false)     ,
	init_hes_cross_done_(false)   ,
	record_newton_atom_done_(false) ,
	init_ran_obj_done_(false)     ,
	init_hes_ran_obj_done_(false)
	{ }
/* $$
$head initialize$$
Directly after construction, use this function to initialize
the derived class object; see $cref/initialize/initialize/$$.
$codep */
	std::map<std::string, size_t> initialize(
		const d_vector& fixed_vec, const d_vector& random_vec
	);
/* $$
$head optimize_random$$
Given the fixed effects, optimize with respect to the random effects;
see  $cref/optimize_random/optimize_random/$$.
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
see  $cref/optimize_fixed/optimize_fixed/$$.
$codep */
	d_vector optimize_fixed(
		const std::string& fixed_options    ,
		const std::string& random_options   ,
		const d_vector&    fixed_lower      ,
		const d_vector&    fixed_upper      ,
		const d_vector&    fix_constraint_lower ,
		const d_vector&    fix_constraint_upper ,
		const d_vector&    fixed_in         ,
		const d_vector&    random_lower     ,
		const d_vector&    random_upper     ,
		const d_vector&    random_in
	);
/* $$
$childtable%src/derived_ctor.omh
	%src/ran_likelihood.omh
	%src/fix_likelihood.omh
	%src/fix_constraint.omh
	%src/initialize.cpp
	%src/optimize_random.cpp
	%src/optimize_fixed.cpp
%$$


$end
*/
private:
/*
------------------------------------------------------------------------------
$begin private$$
$spell
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
$$

$section cppad_mixed: Private Declarations$$
These $code cppad_mixed$$ class declarations are $code private$$.
They are $bold not$$ part of the user API,
they may change with time, and they
can $bold not$$ be used by a derived class object
$cref/mixed_object/derived_ctor/mixed_object/$$.

$childtable%include/cppad/mixed/pack.hpp
	%include/cppad/mixed/unpack.hpp
	%src/init_fix_con.cpp
	%src/init_fix_like.cpp
	%src/init_hes_cross.cpp
	%src/init_hes_ran.cpp
	%src/init_hes_ran_obj.cpp
	%src/init_ran_like.cpp
	%src/init_ran_obj.cpp
	%src/fix_con_eval.cpp
	%src/fix_con_hes.cpp
	%src/fix_con_jac.cpp
	%src/fix_like_eval.cpp
	%src/fix_like_hes.cpp
	%src/fix_like_jac.cpp
	%src/logdet_jac.cpp
	%src/ran_like_jac.cpp
	%src/ran_obj_eval.cpp
	%src/ran_obj_jac.cpp
	%src/ran_obj_hes.cpp
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
when $cref/optimizing fixed effects/optimize_fixed/$$.
$codep */
	const bool quasi_fixed_;
/* $$
$head initialize_done_$$
The following flag is false after construction and true after
the corresponding member function is called:
$codep */
	bool                initialize_done_;
	bool                init_fix_like_done_;
	bool                init_fix_con_done_;
	// only called when n_random_ > 0
	bool                init_ran_like_done_;
	bool                init_hes_ran_done_;
	bool                init_hes_cross_done_;
	// only called when n_random_ > 0 and quasi_fixed_ is false
	bool                record_newton_atom_done_;
	bool                init_ran_obj_done_;
	bool                init_hes_ran_obj_done_;
/* $$

$head ran_likelihood$$
If $icode%n_random_% > 0%$$ and $code init_ran_like_done_$$,
$cref/ran_like_fun_/init_ran_like/ran_like_fun_/$$ and
$cref/init_ran_like_a1fun_/init_ran_like/init_ran_like_a1fun_/$$ are
recordings of the user's $cref ran_likelihood$$.
function.
$codep */
	CppAD::ADFun<double>      ran_like_fun_;
	CppAD::ADFun<a1_double>   init_ran_like_a1fun_;
/* $$
The following objects hold information for computing derivatives
with these ADFun objects:

$subhead hes_ran_$$
If $icode%n_random_% > 0%$$ and $code init_hes_ran_done_$$,
$cref/hes_ran_/init_hes_ran/hes_ran_/$$ contains
information for the Hessian of the
$cref/random likelihood
	/theory
	/Random Likelihood, f(theta, u)
/$$
with respect to the random effects; i.e.
$latex f_{u,u} ( \theta , u )$$.
$codep */
	CppAD::mixed::sparse_hes_info hes_ran_;
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
	/theory
	/Random Likelihood, f(theta, u)
/$$
; i.e.  $latex f_{u,\theta} ( \theta , u )$$.
$codep */
	CppAD::mixed::sparse_hes_info hes_cross_;
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
	CppAD::mixed::newton_step   newton_atom_;
/* $$

$head ran_obj_fun_$$
If $icode%n_random_% > 0%$$, quasi_fixed_ is false, and
$code init_ran_obj_done_$$,
this is a recording of the second approximation for the
random part of the Laplace approximation, $latex H( \beta , \theta , u)$$;
see $cref/ran_obj_fun_/init_ran_obj/ran_obj_fun_/$$.
$codep */
	CppAD::ADFun<double>        ran_obj_fun_;   // for computing H_beta_beta
/* $$
The following objects hold information for computing derivatives
with this ADFun object:

$subhead hes_ran_obj_$$
If $icode%n_random_% > 0%$$, quasi_fixed_ is false, and
$code init_hes_ran_obj_done_$$,
$cref/hes_ran_obj_/init_hes_ran_obj/hes_ran_obj_/$$ contains
information for the Hessian of the
$cref/random objective
	/theory
	/Objective
	/Random Objective, r(theta)
/$$
$codep */
	CppAD::mixed::sparse_hes_info hes_ran_obj_;
/* $$

$head fix_likelihood_$$
$cref/fix_like_fun_/init_fix_like/fix_like_fun_/$$
is a recording of the fixed part of the likelihood function; see,
$cref fix_likelihood$$.
$codep */
	CppAD::ADFun<double>        fix_like_fun_;     // g(theta)
/* $$
The following objects hold information for computing derivatives
with this ADFun object:

$subhead fix_like_jac_$$
$cref/fix_like_jac_/init_fix_like/fix_like_jac_/$$
contains information for the Jacobian of the
$cref/fixed likelihood/theory/Fixed Likelihood, g(theta)/$$.
$codep */
	CppAD::mixed::sparse_jac_info fix_like_jac_;
/* $$

$subhead fix_like_hes_$$
If $icode quasi_fixed$$ is false,
$cref/fix_like_hes_/init_fix_like/fix_like_hes_/$$
contains information for the Hessian of the
$cref/fixed likelihood/theory/Fixed Likelihood, g(theta)/$$.
$codep */
	CppAD::mixed::sparse_hes_info fix_like_hes_;
/* $$

$head fix_con_fun_$$
$cref/fix_con_fun_/init_fix_con/fix_con_fun_/$$
is a recording of the fixed part of the likelihood function; see,
$cref fix_constraint$$.
$codep */
	CppAD::ADFun<double>        fix_con_fun_;     // c(theta)
/* $$
The following objects hold information for computing derivatives
with this ADFun object:

$subhead fix_con_jac_$$
$cref/fix_con_jac_/init_fix_con/fix_con_jac_/$$
contains information for the Jacobian of the
constraint function $latex c ( \theta )$$.
$codep */
	CppAD::mixed::sparse_jac_info fix_con_jac_;
/* $$

$subhead fix_con_hes_$$
If $icode quasi_fixed$$ is false,
$cref/fix_con_hes_/init_fix_con/fix_con_hes_/$$
contains information for the Hessian of the
$cref/constraints/fix_constraint/$$ function $latex c( \theta )$$.
The corresponding ADFun object is
$cref/fix_con_fun_/init_fix_con/fix_con_fun_/$$.
$codep */
	CppAD::mixed::sparse_hes_info fix_con_hes_;
/* $$
$head Template Member Functions$$
$subhead pack$$
See $cref pack$$.
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
$subhead unpack$$
See $cref unpack$$.
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
$head Initialization Member Functions$$

$subhead init_fix_con$$
See $cref init_fix_con$$.
$codep */
	void init_fix_con(const d_vector& fixed_vec);
/* $$

$subhead init_fix_like$$
See $cref init_fix_like$$.
$codep */
	void init_fix_like(const d_vector& fixed_vec);
/* $$

$subhead init_hes_cross$$
See $cref init_hes_cross$$.
$codep */
	void init_hes_cross(
		const d_vector& fixed_vec ,
		const d_vector& random_vec
	);
/* $$

$subhead init_hes_ran$$
See $cref init_hes_ran$$.
$codep */
	void init_hes_ran(
		const d_vector& fixed_vec ,
		const d_vector& random_vec
	);
/* $$

$subhead init_hes_ran_obj$$
See $cref init_hes_ran_obj$$.
$codep */
	void init_hes_ran_obj(
		const d_vector& fixed_vec ,
		const d_vector& random_vec
	);
/* $$

$subhead init_ran_like$$
See $cref init_ran_like$$.
$codep */
	void init_ran_like(
		const d_vector& fixed_vec ,
		const d_vector& random_vec
	);
/* $$

$subhead init_ran_obj$$
See $cref init_ran_obj$$.
$codep */
	void init_ran_obj(
		const d_vector& fixed_vec ,
		const d_vector& random_vec
	);
/* $$
$head Other Member Functions$$

$subhead fix_con_eval$$
See $cref fix_con_eval$$
$codep */
	d_vector fix_con_eval(const d_vector& fixed_vec);
	friend bool ::fix_con_eval_xam(void);
/* $$

$subhead fix_con_jac$$
See $cref fix_con_jac$$
$codep */
	void fix_con_jac(
		const d_vector&        fixed_vec   ,
		CppAD::vector<size_t>& row_out     ,
		CppAD::vector<size_t>& col_out     ,
		d_vector&              val_out
	);
	friend bool ::fix_con_jac_xam(void);
/* $$

$subhead fix_con_hes$$
See $cref fix_con_hes$$
$codep */
	void fix_con_hes(
		const d_vector&        fixed_vec   ,
		const d_vector&        weight      ,
		CppAD::vector<size_t>& row_out     ,
		CppAD::vector<size_t>& col_out     ,
		d_vector&              val_out
	);
	friend bool ::fix_con_hes_xam(void);
/* $$

$subhead fix_like_eval$$
See $cref fix_like_eval$$
$codep */
	d_vector fix_like_eval(const d_vector& fixed_vec);
	friend bool ::fix_like_eval_xam(void);
/* $$

$subhead fix_like_jac$$
See $cref fix_like_jac$$
$codep */
	void fix_like_jac(
		const d_vector&        fixed_vec   ,
		CppAD::vector<size_t>& row_out     ,
		CppAD::vector<size_t>& col_out     ,
		d_vector&              val_out
	);
	friend bool ::fix_like_jac_xam(void);
/* $$

$subhead fix_like_hes$$
See $cref fix_like_hes$$
$codep */
	void fix_like_hes(
		const d_vector&        fixed_vec   ,
		const d_vector&        weight      ,
		CppAD::vector<size_t>& row_out     ,
		CppAD::vector<size_t>& col_out     ,
		d_vector&              val_out
	);
	friend bool ::fix_like_hes_xam(void);
/* $$

$subhead logdet_jac$$
See $cref logdet_jac$$
$codep */
	void logdet_jac(
		const d_vector& fixed_vec  ,
		const d_vector& random_vec ,
		d_vector&       logdet_fix ,
		d_vector&       logdet_ran
	);
	friend bool ::logdet_jac_xam(void);
/* $$

$subhead ran_like_jac$$
See $cref ran_like_jac$$
$codep */
	a1d_vector ran_like_jac(
		const a1d_vector&       fixed_vec   ,
		const a1d_vector&       random_vec
	);
	friend bool ::ran_like_jac_xam(void);
/* $$

$subhead ran_obj_eval$$
See $cref ran_obj_eval$$
$codep */
	double ran_obj_eval(
		const d_vector& fixed_vec  ,
		const d_vector& random_vec
	);
	friend bool ::ran_obj_eval_xam(void);
/* $$

$subhead ran_obj_jac$$
See $cref ran_obj_jac$$
$codep */
	void ran_obj_jac(
		const d_vector& fixed_vec  ,
		const d_vector& random_vec ,
		d_vector&       r_fixed
	);
	friend bool ::ran_obj_jac_xam(void);
	friend bool ::der_var_hes(void);
	friend bool ::delta_ran_obj(void);
/* $$

$subhead ran_obj_hes$$
See $cref ran_obj_hes$$
$codep */
	void ran_obj_hes(
		const d_vector&         fixed_vec   ,
		const d_vector&         random_vec  ,
		CppAD::vector<size_t>&  row_out     ,
		CppAD::vector<size_t>&  col_out     ,
		d_vector&               val_out
	);
	friend bool ::ran_obj_hes_xam(void);
/* $$
$end
-------------------------------------------------------------------------------
*/
};


# include <cppad/mixed/pack.hpp>
# include <cppad/mixed/unpack.hpp>

# endif
