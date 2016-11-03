// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
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
	CppAD
	checkpointing
	algo
	bool
	adfun
	const
	hes
	cholesky
	Hlow
	Eigen
$$

$section Newton Step Algorithm Constructor$$

$head Syntax$$
$codei%CppAD::mixed::newton_step_algo %algo%(
	%a1_adfun%, %hes_info%, %theta%, %u%, %cholesky%
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

$head a1_adfun$$
This is a recording of the function $latex f( \theta , u)$$
for which we are checkpointing the Newton step and log determinant for.
The routine $cref/pack(theta, u)/pack/$$ is used to
convert the pair of vectors into the argument vector for $icode a1_adfun$$.

$head hes_info$$
This argument has prototype
$code%
	const CppAD::mixed::sparse_hes_info& %hes_info%
%$$
It is the $cref sparse_hes_info$$ for the Hessian with respect to the
random effects (as a function of the fixed and random effects); i.e.
$latex f_uu ( \theta , u)$$.
The vector $icode%hes_info%.val%$$ is not used.

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

$head cholesky$$
This argument has prototype
$codei%
	CppAD::mixed::sparse_ad_cholesky& %cholesky%
%$$
The member $code cholesky_$$
is set to be a reference to $icode cholesky$$.
In addition,
$codei%
	%cholesky%.initialize( %a1_Hlow% )
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

$head a1_adfun_$$
This member variable has prototype
$codei%
	CppAD::ADFun<a1_double>&          a1_adfun_
%$$
and is a reference to the $icode a1_adfun$$ argument.

$head hes_info_$$
This member variable has prototype
$codei%
	sparse_hes_info   hes_info_;
%$$
see $cref sparse_hes_info$$.
It is set so that a subsequent call of the form
$codei%
	%a1_adfun_%.SparseHessian(
		%a1_theta_u%,
		%a1_w%,
		%not_used%,
		hes_info_.row,
		hes_info_.col,
		%a1_val_out%,
		hes_info_.work
	)
%$$
will compute the Hessian $latex f_{uu} ( \theta , u )$$.
Here the vectors $icode a1_theta_u$$, $icode a1_w$$, and
$icode a1_val_out$$ have prototype
$codei
	CppAD::vector< CppAD::AD<double> > %a1_theta_u%, %a1_w%, %a1_val_out%
%$$
The vector $icode a1_theta_u$$ specifies $latex ( \theta , u )$$,
$icode a1_w$$ has length one and its element has the value one,
$icode a1_val_out$$ is the value of the Hessian at the row and column
indices corresponding to $code hes_info_$$.

$subhead Remark$$
The vector
$cref/hes_info_.val/sparse_hes_info/Sparse Hessian Call/hes_info.val/$$
is not used because the call to $code SparseHessian$$
uses $code CppAD::AD<double>$$ and not $code double$$ values.

$end
*/
// -------------------------------------------------------------------------
// BEGIN PROTOTYPE
newton_step_algo::newton_step_algo(
	CppAD::ADFun<a1_double>&      a1_adfun      ,
	const sparse_hes_info&        hes_info      ,
	const CppAD::vector<double>&  theta         ,
	const CppAD::vector<double>&  u             ,
	sparse_ad_cholesky&           cholesky      )
// END PROTOTYPE
:
n_fixed_ ( theta.size() ) ,
n_random_( u.size()     ) ,
a1_adfun_( a1_adfun     ) ,
hes_info_( hes_info     ) ,
cholesky_( cholesky     )
{	assert( a1_adfun.Domain() == n_fixed_ + n_random_ );
	assert( a1_adfun.Range()  == 1 );
	assert( hes_info.val.size() == 0 );
	//
# if CPPAD_MIXED_USE_ATOMIC_CHOLESKY
	// total number of variables
	size_t n_both = n_fixed_ + n_random_;
	//
	//	create an a1d_vector containing (theta, u)
	a1d_vector a1_theta_u(n_both);
	for(size_t j = 0; j < n_fixed_; j++)
		a1_theta_u[j] = theta[j];
	for(size_t j = 0; j < n_random_; j++)
		a1_theta_u[n_fixed_ + j] = u[j];
	//
	// sparsity pattern not needed once we have hes_info_.work
	CppAD::vectorBool not_used;
	//
	// compute the sparse Hessian
	a1d_vector a1_w(1), a1_val_out( hes_info_.row.size() );
	a1_w[0] = 1.0;
	a1_adfun_.SparseHessian(
		a1_theta_u,
		a1_w,
		not_used,
		hes_info_.row,
		hes_info_.col,
		a1_val_out,
		hes_info_.work
	);
	// set a1_hessian to lower triangular sparse matrix representation of
	// f_uu (theta, u)
	typedef Eigen::SparseMatrix<a1_double, Eigen::ColMajor> a1_eigen_sparse;
	a1_eigen_sparse a1_hessian(n_random_, n_random_);
	size_t K = hes_info_.row.size();
	for(size_t k = 0; k < K; k++)
	{	assert( n_fixed_ <= hes_info_.col[k]  );
		assert( hes_info_.col[k]  <= hes_info_.row[k] );
		size_t i = hes_info_.row[k] - n_fixed_;
		size_t j = hes_info_.col[k] - n_fixed_;
		a1_hessian.insert(i, j) = a1_val_out[k];
	}
	// initialize cholesky_
	cholesky_.initialize( a1_hessian );
# endif
}
/*
------------------------------------------------------------------------------
$begin newton_step_algo$$
$spell
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
	const a1d_vector& a1_theta_u_v    ,
	a1d_vector&       a1_logdet_step  )
// END PROTOTYPE
{	assert( a1_theta_u_v.size() == n_fixed_ + 2 * n_random_ );
	assert( a1_logdet_step.size() == 1 + n_random_ );

	// extract theta_u subvector
	a1d_vector a1_theta_u(n_fixed_ + n_random_);
	for(size_t j = 0; j < n_fixed_ + n_random_; j++)
		a1_theta_u[j] = a1_theta_u_v[j];

	// sparsity pattern not needed once we have hes_info_.work
	CppAD::vectorBool not_used;

	// compute the sparse Hessian
	a1d_vector a1_w(1), a1_val_out( hes_info_.row.size() );
	a1_w[0] = 1.0;
	a1_adfun_.SparseHessian(
		a1_theta_u,
		a1_w,
		not_used,
		hes_info_.row,
		hes_info_.col,
		a1_val_out,
		hes_info_.work
	);

	// declare eigen matrix types
	typedef Eigen::Matrix<a1_double, Eigen::Dynamic, 1>     a1_eigen_vector;
	typedef Eigen::SparseMatrix<a1_double, Eigen::ColMajor> a1_eigen_sparse;

	// set a1_hessian to lower triangular eigen sparse matrix representation of
	// f_uu (theta, u)
	a1_eigen_sparse a1_hessian(n_random_, n_random_);
	size_t K = hes_info_.row.size();
	for(size_t k = 0; k < K; k++)
	{	assert( n_fixed_ <= hes_info_.col[k]  );
		assert( hes_info_.col[k]  <= hes_info_.row[k] );
		size_t i = hes_info_.row[k] - n_fixed_;
		size_t j = hes_info_.col[k] - n_fixed_;
		a1_hessian.insert(i, j) = a1_val_out[k];
	}

	// set rhs to an eigen vector representation of v
	a1_eigen_vector rhs(n_random_);
	for(size_t j = 0; j < n_random_; j++)
		rhs(j) = a1_theta_u_v[n_fixed_ + n_random_ + j];

	// compute the newt step and the log determinant
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
newton_step::newton_step(void)
: checkpoint_fun_(CPPAD_MIXED_NULL_PTR)
{ }
// newton_step destructor
newton_step::~newton_step(void)
{	if( checkpoint_fun_ != CPPAD_MIXED_NULL_PTR )
		delete checkpoint_fun_;
}
/*
-------------------------------------------------------------------------------
$begin newton_step_initialize$$
$spell
	bool
	adfun
	CppAD
	checkpoint
	checkpointing
	hes
	const
$$

$section Initialize the Newton Step Checkpoint Function$$

$head Syntax$$
$icode%newton_checkpoint%.initialize(%a1_adfun%, %hes_info%, %theta%, %u%)
%$$

$head Prototype$$
$srcfile%src/eigen/newton_step.cpp
	%4%// BEGIN PROTOTYPE%// END PROTOTYPE%5%$$

$head Private$$
This class is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head newton_checkpoint$$
This is a Newton step object; see
$cref/newton_checkpoint/newton_step_ctor/newton_checkpoint/$$.

$head a1_adfun$$
This is a recording of the function $latex f( \theta , u)$$
for which we are checkpointing the Newton step and log determinant for.
The routine $cref/pack(theta, u)/pack/$$ is used to
convert the pair of vectors into the argument vector for $icode a1_adfun$$.

$head hes_info$$
This argument has prototype
$code%
	const CppAD::mixed::sparse_hes_info& %hes_info%
%$$
It is the $cref sparse_hes_info$$ for the Hessian with respect to the
random effects (as a function of the fixed and random effects); i.e.
$latex f_uu ( \theta , u)$$.
The vector $icode%hes_info%.val%$$ is not used.

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
	CppAD::ADFun<a1_double>&          a1_adfun      ,
	const sparse_hes_info&            hes_info      ,
	const CppAD::vector<double>&      theta         ,
	const CppAD::vector<double>&      u             )
// END PROTOTYPE
{	assert( checkpoint_fun_ == CPPAD_MIXED_NULL_PTR );
	//
	size_t n_fixed  = theta.size();
	size_t n_random = u.size();
	// algo
	newton_step_algo algo( a1_adfun, hes_info, theta, u, cholesky_);
	// checkpoint_fun_
	a1d_vector a1_theta_u_v(n_fixed + 2 * n_random);
	for(size_t j = 0; j < n_fixed; j++)
		a1_theta_u_v[j] = theta[j];
	for(size_t j = 0; j < n_random; j++)
	{	a1_theta_u_v[n_fixed + j] = u[j];
		// right hand size should not matter
		a1_theta_u_v[n_fixed + n_random + j] = 1.0;
	}
	const char* name = "cppad_mixed newton_step";
	a1d_vector a1_logdet_step(1 + n_random);
	CppAD::atomic_base<double>::option_enum sparsity =
		CppAD::atomic_base<double>::pack_sparsity_enum;
	bool optimize = bool( CPPAD_MIXED_OPTIMIZE_CPPAD_FUNCTION );
	checkpoint_fun_ = new CppAD::checkpoint<double>(
		name, algo, a1_theta_u_v, a1_logdet_step, sparsity, optimize
	);
	assert( checkpoint_fun_ != CPPAD_MIXED_NULL_PTR );
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
	const a1d_vector& a1_theta_u_v  ,
	a1d_vector&       a1_logdet_step )
// END PROTOTYPE
{	assert( checkpoint_fun_ != CPPAD_MIXED_NULL_PTR );
	(*checkpoint_fun_)(a1_theta_u_v, a1_logdet_step);
}

} } // END_CPPAD_MIXED_NAMESPACE
