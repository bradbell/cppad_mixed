// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_NEWTON_STEP_HPP
# define CPPAD_MIXED_NEWTON_STEP_HPP
/*
$begin newton_step$$
$spell
	CppAD
	cppad
	logdet
	adfun
	sv
	var
	eval
	CppAD
	checkpointing
	const
	bool
$$

$section Checkpoint Newton Step and Log Determinant Calculation$$

$head Syntax$$
$codei%CppAD::mixed::newton_step %newton_checkpoint%()
%$$
$icode%newton_checkpoint%.initialize(%bool_sparsity%, %a1_adfun%, %theta%, %u%)
%$$
$icode%sv% = newton_checkpoint%.size_var()
%$$
$icode%newton_checkpoint%.eval(%a1_theta_u_v%, %a1_logdet_step%)
%$$

$head Private$$
This class is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
This stores the operation sequence corresponding to one Newton step
and a log-determinant calculation.

$head Constructor$$
The sparse Hessian checkpoint object constructor
$codei%
	newton_step %newton_checkpoint%()
%$$
creates a $code CppAD::checkpoint<double>$$ object for evaluating
the log of the determinant
$latex \[
( \theta , u ) \rightarrow \log \det [ f_{uu} ( \theta , u )  ]
\] $$
and the Newton step
$latex \[
( \theta , u , v ) \rightarrow f_{uu} ( \theta , u )^{-1} v
\] $$

$head Destructor$$
The object $icode newton_checkpoint$$ must still exist (not be destructed)
for as long as any $code CppAD::ADFun$$ objects use its checkpoint operation.

$head initialize$$
The $icode newton_checkpoint$$ object must be initialized,
before any calls to its $code eval$$ routine, using the syntax
$codei%
	newton_checkpoint.initialize(%bool_sparsity%, %a1_adfun%, %theta%, %u%)
%$$

$head bool_sparsity$$
This argument has prototype
$codei%
	bool %bool_sparsity%
%$$
If this is true, use boolean patterns
(otherwise use set sparsity patterns)
when computing the sparsity
for the Hessian w.r.t the random effects of the
$cref/random likelihood/theory/Random Likelihood, f(theta, u)/$$;
i.e. $latex f_{uu} ( \theta , u )$$.

$head a1_adfun$$
This $code initialize$$ argument has prototype
$codei%
	CppAD::ADFun< CppAD::AD<double> > %a1_adfun%
%$$
This is a recording of the function $latex f( \theta , u)$$
for which we are checkpointing the Newton step and log determinant for.
The routine $cref/pack(theta, u)/pack/$$ is used to
convert the pair of vectors into the argument vector for $icode a1_adfun$$.

$head theta$$
This $code initialize$$ argument has prototype
$codei%
	const CppAD::vector<double>& %theta%
%$$
It is a value for $latex \theta$$
at which we can evaluate the Newton step and
log determinant.

$head u$$
This $code initialize$$ argument has prototype
$codei%
	const CppAD::vector<double>& %u%
%$$
It is a value for $latex u$$
at which we can evaluate the Newton step and
log determinant.

$head size_var$$
The return value for this member function has prototype
$codei%
	size_t %sv%
%$$
It is the number of variables in the tape used to represent the
CppAD checkpoint function.
If the  $code initialize$$ member function has not yet been called,
the return value is $icode%sv% = 0%$$.

$head eval$$
The $code initialize$$ member function must be called before
the $code eval$$ member function is used.

$subhead a1_theta_u_v$$
This $code eval$$ argument has prototype
$codei%
	const CppAD::vector< CppAD::AD<double> >& %a1_theta_u_v%
%$$
It is a point at which we are evaluating the Newton step
and log determinant.
This vector has size $icode%theta%.size() + 2 * %u%.size()%$$
and the order of its elements are $latex \theta$$, followed by $latex u$$
followed by $latex v$$.

$subhead a1_logdet_step$$
This $code eval$$ argument has prototype
$codei%
	CppAD::vector< CppAD::AD<double> >& %a1_logdet_step%
%$$
It size is $codei%1 + %u%.size()%$$.
The input value of its elements does not matter.
Upon return,
$codei%
	%a1_logdet_step%[0] =%$$ $latex \log \det [ f_{uu} ( \theta , u )  ]$$
$pre
$$
and for $icode%j% = 1 ,%...%, %m%%$$,
$codei%
	%a1_logdet_step%[j] =%$$ $latex s_{j-1}$$
$pre
$$
where $latex s$$ is the Newton step; i.e.,
$latex \[
	s = f_{uu} ( \theta , u )^{-1} v
\] $$

$head Example$$
The file $cref newton_step.cpp$$ is an example
and test of $code newton_step$$.

$childtable%example/private/newton_step.cpp
	%src/eigen/newton_step.cpp
%$$

$end
------------------------------------------------------------------------------
*/

# include <cppad/cppad.hpp>
# include <cppad/mixed/sparse_hes_info.hpp>
# include <cppad/mixed/sparse_ad_cholesky.hpp>

// This value 1 (true) is not yet working
# define CPPAD_MIXED_USE_SPARSE_AD_CHOLESKY 0

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE


// class used to evalute algorithm that is checkpointed
class newton_step_algo {
	typedef CppAD::AD<double>            a1_double;
	typedef CppAD::vector<a1_double>     a1d_vector;
private:
		// number of fixed effects
		const size_t                      n_fixed_;
		// niumber of random effects
		const size_t                      n_random_;
		// random likelihood f(theta, u)
		CppAD::ADFun<a1_double>&          a1_adfun_;
		// information for computing f_uu (theta , u)
		// (hes_info.val is not used and has size zero)
		sparse_hes_info                   hes_info_;
		// reference to sparse Cholesky factorization in newton step object
		sparse_ad_cholesky&               cholesky_;
public:
	// constructor for algorithm that is checkpointed
	newton_step_algo(
		bool                              bool_sparsity ,
		CppAD::ADFun<a1_double>&          a1_adfun      ,
		const CppAD::vector<double>&      theta         ,
		const CppAD::vector<double>&      u             ,
		sparse_ad_cholesky&               cholesky
	);
	// evaluates algorithm that is checkpointed
	void operator()(
		const a1d_vector& a1_theta_u_v   ,
		a1d_vector&       a1_logdet_step
	);
};

// class used to hold the checkpoint version of the algorithm
class newton_step {
	typedef CppAD::AD<double>             a1_double;
	typedef CppAD::vector<a1_double>      a1d_vector;
private:
	// checkpoint version of the newton step algorithm
	CppAD::checkpoint<double>*            checkpoint_fun_;
	//
	// Sparse Cholesky factorization as an atomic AD operation.
	// Stored here becasue newton_step_algo object does not persist
	// for the lenght of time that cholesky_ call backs are used by CppAD.
	sparse_ad_cholesky                    cholesky_;
public:
	// constuctor
	newton_step(void);
	// destructor
	~newton_step(void);
	// setup the checkpoint function
	void initialize(
		bool                              bool_sparsity ,
		CppAD::ADFun<a1_double>&          a1_adfun      ,
		const CppAD::vector<double>&      fixed_vec     ,
		const CppAD::vector<double>&      random_vec
	);
	// size of the checkpoint function
	size_t size_var(void);
	// use the checkpoint function
	void eval(const a1d_vector& a1_theta_u_v, a1d_vector& a1_logdet_step);
};


} } // END_CPPAD_MIXED_NAMESPACE

# endif
