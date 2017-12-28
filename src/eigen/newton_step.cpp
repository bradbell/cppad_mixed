// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-17 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <Eigen/Sparse>
# include <cppad/mixed/newton_step.hpp>
# include <cppad/mixed/configure.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

/*
$begin newton_step_algo_ctor$$
$spell
	rcv
	CppAD
	checkpointing
	algo
	bool
	adfun
	const
	hes
	cholesky cholesky
	Hlow
	Eigen
	jac
$$

$section Newton Step Algorithm Constructor$$

$head Syntax$$
$codei%CppAD::mixed::newton_step_algo %algo%(
	%a1fun%, %jac_a1fun%, %hes_rcv%, %hes_work%, %theta%, %u%,
)%$$

$head Prototype$$
$srcfile%src/eigen/newton_step.cpp
	%4%// BEGIN PROTOTYPE%// END PROTOTYPE%1%$$

$head Private$$
This class is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head algo$$
This is the Newton step algorithm object that is constructed
by this operation.

$head a1fun$$
This is a recording of the function $latex f( \theta , u)$$
for which we are checkpointing the Newton step and log determinant for.
The routine $cref/pack(theta, u)/pack/$$ is used to
convert the pair of vectors into the argument vector for $icode a1fun$$.
A reference to $icode a1fun$$ is stored,
so this object should exist as long as $icode algo$$ exists.

$head jac_a1fun$$
This is a recording of the Partial $latex f_u ( \theta , u)$$.
A reference to $icode jac_a1fun$$ is stored,
so this object should exist as long as $icode algo$$ exists.

$head hes_rcv$$
This argument has prototype
$code%
	const CppAD::mixed::sparse_rcv& %hes_rcv%
%$$
It is the
$cref/sparse_rcv/typedef/Sparse Types/sparse_rcv/$$
information for the lower triangle of the Hessian with respect to the
random effects, as a function of the fixed and random effects; i.e.
$latex f_uu ( \theta , u)$$.
In addition, this matrix is in column major order.

$head hes_work$$
This argument has prototype
$code%
	const CppAD::sparse_hes_work& %hes_work%
%$$
It contains the necessary information so that the call
$codei%
	a1fun.sparse_hes(
		%a1_x%,
		%a1_w%
		a1_hes_rcv,
		%not_used_pattern%,
		%not_used_coloring%,
		hes_work
	)
%$$
can be used the compute the sparse Hessian.
Here $icode a1_hes_rcv$$ is the same as $icode hes_rcv$$, except
that its value vector has type $code a1_vector$$ instead of $code d_vector$$.

$head theta$$
This is a value for $latex \theta$$
at which we can evaluate the Newton step and log determinant.
If $cref/CPPAD_MIXED_USE_ATOMIC_CHOLESKY
	/configure.hpp
	/CPPAD_MIXED_USE_ATOMIC_CHOLESKY
/$$
is $code 0$$, $icode theta$$ is not used.

$head u$$
This is a value for $latex u$$
at which we can evaluate the Newton step and log determinant.
If $cref/CPPAD_MIXED_USE_ATOMIC_CHOLESKY
	/configure.hpp
	/CPPAD_MIXED_USE_ATOMIC_CHOLESKY
/$$
is $code 0$$, $icode u$$ is not used.

$head cholesky_$$
This member variable has prototype
$codei%
	CppAD::mixed::sparse_ad_cholesky cholesky_
%$$
$codei%
	cholesky_.initialize( %a1_Hlow% )
%$$
is called where $icode a1_Hlow$$ has prototype
$codei%
	Eigen::SparseMatrix< CppAD::AD<double>, Eigen::ColMajor>& %a1_Hlow%
%$$
and is the lower triangle of the Hessian $latex f_{uu} ( \theta , u )$$.

$head n_fixed_$$
This member variable has prototype
$codei%
	const size_t %n_fixed_%
%$$
and is the number of fixed effects; i.e., the size of $icode theta$$.

$head n_random_$$
This member variable has prototype
$codei%
	const size_t n_random_
%$$
and is the number of random effects; i.e., the size of $icode u$$.

$head a1fun_$$
This member variable has prototype
$codei%
	CppAD::ADFun<a1_double>&          a1fun_
%$$
and is a reference to the $icode a1fun$$ argument.

$head a1_jac2hes_rcv_$$
This member variable has prototype
$codei%
	a1_sparse_rcv_&   a1_jac2hes_rcv_
%$$

$end
*/
// -------------------------------------------------------------------------
// BEGIN PROTOTYPE
newton_step_algo::newton_step_algo(
	CppAD::ADFun<a1_double>&      a1fun         ,
	CppAD::ADFun<a1_double>&      jac_a1fun     ,
	const sparse_rcv&             hes_rcv       ,
	const CppAD::sparse_hes_work& hes_work      ,
	const CppAD::vector<double>&  theta         ,
	const CppAD::vector<double>&  u             )
// END PROTOTYPE
:
n_fixed_  ( theta.size() ) ,
n_random_ ( u.size()     ) ,
a1fun_    ( a1fun        ) ,
jac_a1fun_( jac_a1fun    ) ,
hes_work_ ( hes_work     )
{	assert( a1fun_.Range() == 1 );
	assert( a1fun.Domain() == n_fixed_ + n_random_ );
	//
	// total number of variables
	size_t n_both = n_fixed_ + n_random_;
	//
	//	create an a1_vector containing (theta, u)
	a1_vector a1_theta_u(n_both);
	for(size_t j = 0; j < n_fixed_; j++)
		a1_theta_u[j] = theta[j];
	for(size_t j = 0; j < n_random_; j++)
		a1_theta_u[n_fixed_ + j] = u[j];
	//
	// set a1_jac2hes_rcv_
	sparse_rc pattern(n_random_, n_fixed_ + n_random_, hes_rcv.nnz());
	for(size_t k = 0; k < hes_rcv.nnz(); k++)
	{	size_t r = hes_rcv.row()[k];
		size_t c = hes_rcv.col()[k];
		assert( n_fixed_ <= r );
		assert( n_fixed_ <= c );
		pattern.set(k, r - n_fixed_, c);
	}
	a1_jac2hes_rcv_ = a1_sparse_rcv( pattern );
	for(size_t k = 0; k < hes_rcv.nnz(); k++)
		a1_jac2hes_rcv_.set(k, hes_rcv.val()[k] );
	// =======================================================================
# if CPPAD_MIXED_USE_ATOMIC_CHOLESKY
	// set a1_hessian to lower triangular sparse matrix representation of
	// f_uu (theta, u)
	typedef Eigen::SparseMatrix<a1_double, Eigen::ColMajor> a1_eigen_sparse;
	a1_eigen_sparse a1_hessian(n_random_, n_random_);
	for(size_t k = 0; k < a1_jac2hes_rcv_.nnz(); k++)
	{	size_t r = a1_jac2hes_rcv_.row()[k];
		size_t c = a1_jac2hes_rcv_.col()[k];
		//
		assert( r < n_random_ );
		assert( n_fixed_ <= c );
		assert( c <= r + n_fixed_ );
		//
		size_t i = r;
		size_t j = c - n_fixed_;
		a1_hessian.insert(i, j) = a1_jac2hes_rcv_.val()[k];
	}
	// initialize cholesky_
	cholesky_.initialize( a1_hessian );
# endif
}
/*
------------------------------------------------------------------------------
$begin newton_step_algo$$
$spell
	rcv
	algo
	logdet
	CppAD
$$

$section Newton Step Algorithm Evaluation$$

$head Syntax$$
$icode%algo%(%a1_theta_u_v%, %a1_logdet_step%)%$$

$head Prototype$$
$srcfile%src/eigen/newton_step.cpp
	%4%// BEGIN PROTOTYPE%// END PROTOTYPE%3%$$

$head Private$$
This class is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head algo$$
This is a Newton step object; see
$cref/algo/newton_step_algo_ctor/algo/$$.

$head a1_theta_u_v$$
This is a point at which we are evaluating the Newton step
and log determinant.
This vector has size $code n_fixed_ + 2 * n_random_$$.
and the order of its elements are $latex \theta$$, followed by $latex u$$
followed by $latex v$$.

$head a1_logdet_step$$
This vector has  size $code 1 + n_random_$$.
The input value of its elements does not matter.
Upon return,
$codei%
	%a1_logdet_step%[0]
%$$
is the log of the determinant of $latex f_{uu} ( \theta , u )$$.
For $icode%j% = 1 ,%...%, n_random_%$$,
$codei%
	%a1_logdet_step%[%j%]
%$$
is the $latex s_{j-1}$$ where $latex s$$ is the Newton step; i.e.,
$latex \[
	s = f_{uu} ( \theta , u )^{-1} v
\] $$

$end
*/
// BEGIN PROTOTYPE
void newton_step_algo::operator()(
	const a1_vector&  a1_theta_u_v    ,
	a1_vector&        a1_logdet_step  )
// END PROTOTYPE
{	assert( a1_theta_u_v.size() == n_fixed_ + 2 * n_random_ );
	assert( a1_logdet_step.size() == 1 + n_random_ );

	// extract theta_u subvector
	a1_vector a1_theta_u(n_fixed_ + n_random_);
	for(size_t j = 0; j < n_fixed_ + n_random_; j++)
		a1_theta_u[j] = a1_theta_u_v[j];

	// sparsity pattern not needed once we have hes_work_
	CppAD::vectorBool not_used;

	// compute the sparse Hessian
	jac_a1fun_.subgraph_jac_rev(a1_theta_u, a1_jac2hes_rcv_);

	// declare eigen matrix types
	typedef Eigen::Matrix<a1_double, Eigen::Dynamic, 1>     a1_eigen_vector;
	typedef Eigen::SparseMatrix<a1_double, Eigen::ColMajor> a1_eigen_sparse;

	// set a1_hessian to lower triangular eigen sparse matrix representation of
	// f_uu (theta, u)
	a1_eigen_sparse a1_hessian(n_random_, n_random_);
	for(size_t k = 0; k < a1_jac2hes_rcv_.nnz(); k++)
	{	size_t r = a1_jac2hes_rcv_.row()[k];
		size_t c = a1_jac2hes_rcv_.col()[k];
		//
		assert( r < n_random_ );
		assert( n_fixed_ <= c );
		assert( c <= r + n_fixed_ );
		//
		size_t i = r;
		size_t j = c - n_fixed_;
		a1_hessian.insert(i, j) = a1_jac2hes_rcv_.val()[k];
	}

	// set rhs to an eigen vector representation of v
	a1_eigen_vector rhs(n_random_);
	for(size_t j = 0; j < n_random_; j++)
		rhs(j) = a1_theta_u_v[n_fixed_ + n_random_ + j];

	// compute the newton step and the log determinant
	a1_eigen_vector step;
	a1_double       logdet = 0.0;
# if CPPAD_MIXED_USE_ATOMIC_CHOLESKY
	// Permutation matrix corresponding to this sparse_ad_cholesky
	const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>& P =
		cholesky_.permutation();
	//
	// compute P^T * L * L^T * P = f_uu (theta, u) factorization
	a1_eigen_sparse a1_L;
	cholesky_.eval(a1_hessian, a1_L);
	//
	// compute logdet of f_uu (theta, u)
	for(size_t j = 0; j < n_random_; j++)
	{	a1_eigen_sparse::InnerIterator itr(a1_L, j);
		// first entry in each column is the diagonal entry
		assert( itr.row() == itr.col() );
		logdet += log( itr.value() );
	}
	logdet *= 2.0;
	//
	// solve the equation f_uu (theta , u) * step = v
	// P^T * L * L^T * P * step = v
	a1_eigen_sparse a1_U = a1_L.transpose();
	a1_eigen_vector tmp  = P * rhs;
	tmp  = a1_L.triangularView<Eigen::Lower>().solve(tmp);
	tmp  = a1_U.triangularView<Eigen::Upper>().solve(tmp);
	step = P.transpose() * tmp;
# else
	// compute an L * D * L^T Cholesky factorization of f_{u,u}(theta, u)
	Eigen::SimplicialLDLT<a1_eigen_sparse, Eigen::Lower> chol;
	chol.analyzePattern(a1_hessian);
	chol.factorize(a1_hessian);

	// solve the equation f_{u,u} ( theta, u ) * step = v
	step = chol.solve(rhs);

	// compute the logdet( f_{u,u} (theta, u ) )
	a1_eigen_vector diag = chol.vectorD();
	for(size_t j = 0; j < n_random_; j++)
		logdet += log( diag(j) );
# endif

	// pack resutls into return vector
	a1_logdet_step[0] = logdet;
	for(size_t j = 0; j < n_random_; j++)
		a1_logdet_step[1 + j] = step(j);

	return;
}
/*
-------------------------------------------------------------------------------
$begin newton_step_ctor$$
$spell
	CppAD
	eval
	algo
$$

$section Newton Step Checkpoint Function Constructor$$

$head Syntax$$
$codei%CppAD::mixed::newton_step %newton_checkpoint%()
%$$

$head Private$$
This class is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head newton_checkpoint$$
This is the CppAD checkpoint function object that is constructed
by this operation.

$head algo_$$
This member variable has prototype
$codei%
	CppAD::mixed::newton_step_algo* algo_;
%$$
It is set to $code null$$ by the constructor.
It must constructed during the
$cref/initialize/newton_step/initialize/$$ process.

$head checkpoint_fun_$$
This member variable has prototype
$codei%
	CppAD::checkpoint<double>*  checkpoint_fun_;
%$$
It is set to $code null$$ by the constructor.
This ensures that
$cref/initialize/newton_step/initialize/$$
is called before
$cref/eval/newton_step/eval/$$ is used.

$end
*/
// newton_step ctor
newton_step::newton_step(void) :
algo_          (CPPAD_MIXED_NULL_PTR)   ,
checkpoint_fun_(CPPAD_MIXED_NULL_PTR)
{ }
// newton_step destructor
newton_step::~newton_step(void)
{	if( algo_ != CPPAD_MIXED_NULL_PTR )
		delete algo_;
	if( checkpoint_fun_ != CPPAD_MIXED_NULL_PTR )
		delete checkpoint_fun_;
}
/*
-------------------------------------------------------------------------------
$begin newton_step_initialize$$
$spell
	rcv
	bool
	adfun
	CppAD
	checkpoint
	checkpointing
	hes
	const
	jac
$$

$section Initialize the Newton Step Checkpoint Function$$

$head Syntax$$
$icode%newton_checkpoint%.initialize(
	%a1fun%, %hes_rcv%, %hes_work%, %theta%, %u%
)%$$

$head Prototype$$
$srcfile%src/eigen/newton_step.cpp
	%4%// BEGIN PROTOTYPE%// END PROTOTYPE%5%$$

$head Private$$
This class is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head newton_checkpoint$$
This is a Newton step object; see
$cref/newton_checkpoint/newton_step_ctor/newton_checkpoint/$$.

$head a1fun$$
This is a recording of the function $latex f( \theta , u)$$
for which we are checkpointing the Newton step and log determinant for.
The routine $cref/pack(theta, u)/pack/$$ is used to
convert the pair of vectors into the argument vector for $icode a1fun$$.
A reference to $icode a1fun$$ is stored,
so this object should exist as long as $icode newton_checkpoint$$ exists.

$head jac_a1fun$$
This is a recording of the Partial $latex f_u ( \theta , u)$$.
A reference to $icode jac_a1fun$$ is stored,
so this object should exist as long as $icode newton_checkpoint$$ exists.

$head hes_rcv$$
This argument has prototype
$code%
	const CppAD::mixed::sparse_rcv& %hes_rcv%
%$$
It is the
$cref/sparse_rcv/typedef/Sparse Types/sparse_rcv/$$
information for the Hessian with respect to the
random effects (as a function of the fixed and random effects); i.e.
$latex f_uu ( \theta , u)$$.

$head hes_work$$
This argument has prototype
$code%
	const CppAD::sparse_hes_work& %hes_work%
%$$
It contains the necessary information so that the call
$codei%
	a1fun.sparse_hes(
		%a1_x%,
		%a1_w%
		a1_hes_rcv,
		%not_used_pattern%,
		%not_used_coloring%,
		hes_work
	)
%$$
can be used the compute the sparse Hessian.
Here $icode a1_hes_rcv$$ is the same as $icode hes_rcv$$, except
that its value vector has type $code a1_vector$$ instead of $code d_vector$$.

$head theta$$
This is a value for $latex \theta$$
at which we can evaluate the Newton step and log determinant.

$head u$$
This is a value for $latex u$$
at which we can evaluate the Newton step and log determinant.

$head checkpoint_fun_$$
This member variable has prototype
$codei%
	CppAD::checkpoint<double>*  checkpoint_fun_;
%$$
It must be $code null$$ when this function is called.
Upon return,
it will point to a CppAD checkpoint function that is used
to evaluate the Newton step; see $cref newton_step_eval$$.


$end
*/
// BEGIN PROTOTYPE
void newton_step::initialize(
	CppAD::ADFun<a1_double>&          a1fun         ,
	CppAD::ADFun<a1_double>&          jac_a1fun     ,
	const sparse_rcv&                 hes_rcv       ,
	const CppAD::sparse_hes_work      hes_work      ,
	const CppAD::vector<double>&      theta         ,
	const CppAD::vector<double>&      u             )
// END PROTOTYPE
{	assert( algo_           == CPPAD_MIXED_NULL_PTR );
	assert( checkpoint_fun_ == CPPAD_MIXED_NULL_PTR );
	// create algo
	//
	algo_ = new newton_step_algo(
		a1fun, jac_a1fun, hes_rcv, hes_work, theta, u
	);
	assert( algo_ != CPPAD_MIXED_NULL_PTR );
	//
# if CPPAD_MIXED_CHECKPOINT_NEWTON_STEP
	//  create checkpoint_fun_
	//
	size_t n_fixed  = theta.size();
	size_t n_random = u.size();
	a1_vector a1_theta_u_v(n_fixed + 2 * n_random);
	for(size_t j = 0; j < n_fixed; j++)
		a1_theta_u_v[j] = theta[j];
	for(size_t j = 0; j < n_random; j++)
	{	a1_theta_u_v[n_fixed + j] = u[j];
		// right hand size should not matter
		a1_theta_u_v[n_fixed + n_random + j] = 1.0;
	}
	const char* name = "cppad_mixed newton_step";
	a1_vector a1_logdet_step(1 + n_random);
	CppAD::atomic_base<double>::option_enum sparsity =
		CppAD::atomic_base<double>::pack_sparsity_enum;
	bool optimize = bool( CPPAD_MIXED_OPTIMIZE_CPPAD_FUNCTION );
	checkpoint_fun_ = new CppAD::checkpoint<double>(
		name, *algo_, a1_theta_u_v, a1_logdet_step, sparsity, optimize
	);
	assert( checkpoint_fun_ != CPPAD_MIXED_NULL_PTR );
# endif
}
/*
-----------------------------------------------------------------------------
$begin newton_step_size_var$$
$spell
	CppAD
$$

$section Number of Variables in Checkpoint Function$$

$head Private$$
This class is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$srccode%cpp% */
size_t newton_step::size_var(void)
{	if( checkpoint_fun_ ==  CPPAD_MIXED_NULL_PTR )
		return 0;
	return checkpoint_fun_->size_var();
}
/* %$$
$end
-------------------------------------------------------------------------------
$begin newton_step_eval$$
$spell
	eval
	logdet
	CppAD
$$

$section Using the Newton Step Checkpoint Function$$

$head Syntax$$
$icode%newton_checkpoint%.eval(%a1_theta_u_v%, %a1_logdet_step%)
%$$

$head Prototype$$
$srcfile%src/eigen/newton_step.cpp
	%4%// BEGIN PROTOTYPE%// END PROTOTYPE%7%$$

$head Private$$
This class is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Restriction$$
The $cref newton_step_initialize$$ routine must be called
before this routine be used.

$head newton_checkpoint$$
This is a Newton step object; see
$cref/newton_checkpoint/newton_step_ctor/newton_checkpoint/$$.

$head a1_theta_u_v$$
This is a point at which we are evaluating the Newton step
and log determinant.
This vector has size $icode%n_fixed% + 2 * %n_random%$$
and the order of its elements are $latex \theta$$, followed by $latex u$$
followed by $latex v$$.

$head a1_logdet_step$$
This vector has  size $codei%1 + %n_random%$$.
The input value of its elements does not matter.
Upon return,
$codei%
	%a1_logdet_step%[0]
%$$
is the log of the determinant of $latex f_{uu} ( \theta , u )$$.
For $icode%j% = 1 ,%...%, n_random_%$$,
$codei%
	%a1_logdet_step%[%j%]
%$$
is the $latex s_{j-1}$$ where $latex s$$ is the Newton step; i.e.,
$latex \[
	s = f_{uu} ( \theta , u )^{-1} v
\] $$

$end
*/
// BEGIN PROTOTYPE
void newton_step::eval(
	const a1_vector& a1_theta_u_v  ,
	a1_vector&        a1_logdet_step )
// END PROTOTYPE
{	assert( algo_ != CPPAD_MIXED_NULL_PTR );
	//
# if CPPAD_MIXED_CHECKPOINT_NEWTON_STEP
	assert( checkpoint_fun_ != CPPAD_MIXED_NULL_PTR );
	(*checkpoint_fun_)(a1_theta_u_v, a1_logdet_step);
# else
	(*algo_)(a1_theta_u_v, a1_logdet_step);
# endif
}

} } // END_CPPAD_MIXED_NAMESPACE
