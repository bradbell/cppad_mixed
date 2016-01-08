// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
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
$$

$section Checkpoint Newton Step and Log Determinant Calculation$$

$head Syntax$$
$codei%CppAD::mixed::newton_step %newton_atom%()
%$$
$icode%newton_atom%.initialize(%a1_adfun%, %theta%, %u%)
%$$
$icode%sv% = newton_atom%.size_var()
%$$
$icode%newton_atom%.eval(%a1_theta_u_v%, %a1_logdet_step%)
%$$

$head Purpose$$
This stores the operation sequence corresponding to the mappings

$head Constructor$$
The sparse Hessian checkpoint object constructor
$codei%
	newton_step %newton_atom%()
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
The object $icode newton_atom$$ must still exist (not be destructed)
for as long as any $code CppAD::ADFun$$ objects use its atomic operation.

$head initialize$$
The $icode newton_atom$$ object must be initialized,
before any calls to its $code eval$$ routine, using the syntax
$codei%
	newton_atom.initialize(%a1_adfun%, %theta%, %u%)
%$$

$subhead a1_adfun$$
This $code initialize$$ argument has prototype
$codei%
	CppAD::ADFun< CppAD::AD<double> > %a1_adfun%
%$$
This is a recording of the function $latex f(theta, u)$$
for which we are checkpointing the Newton step and log determinant for.
The routine $cref/pack(theta, u)/pack/$$ is used to
convert the pair of vectors into the argument vector for $icode a1_adfun$$.

$subhead theta$$
This $code initialize$$ argument has prototype
$codei%
	const CppAD::vector<double>& %theta%
%$$
It is a value for $latex \theta$$
at which we can evaluate the Newton step and
log determinant.

$subhead u$$
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

$children%example/private/newton_step_xam.cpp
%$$
$head Example$$
The file $cref newton_step_xam.cpp$$ is an example
and test of $code newton_step$$.


$end
----------------------------------------------------------------------------
*/
# include <Eigen/Sparse>
# include <cppad/mixed/newton_step.hpp>
# include <cppad/mixed/configure.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// -------------------------------------------------------------------------
// newton_step_alog ctor
newton_step_algo::newton_step_algo(
	CppAD::ADFun<a1_double>&      a1_adfun ,
	const CppAD::vector<double>&  theta    ,
	const CppAD::vector<double>&  u        )
:
n_fixed_ ( theta.size() ) ,
n_random_( u.size()     ) ,
a1_adfun_( a1_adfun     )
{	assert( row_.size() == 0 );
	assert( col_.size() == 0 );
	assert( a1_adfun.Domain() == n_fixed_ + n_random_ );
	assert( a1_adfun.Range()  == 1 );
	//
	// total number of variables
	size_t n_both = n_fixed_ + n_random_;

	//	create an a1d_vector containing (theta, u)
	a1d_vector a1_theta_u(n_both);
	for(size_t j = 0; j < n_fixed_; j++)
		a1_theta_u[j] = theta[j];
	for(size_t j = 0; j < n_random_; j++)
		a1_theta_u[n_fixed_ + j] = u[j];

	// compute Jacobian sparsity corresponding to parital w.r.t. random effects
	typedef CppAD::vectorBool sparsity_pattern;
	sparsity_pattern r(n_both * n_both);
	for(size_t i = 0; i < n_both; i++)
	{	for(size_t j = 0; j < n_both; j++)
			r[i * n_both + j] = (i >= n_fixed_) && (i == j);
	}
	a1_adfun_.ForSparseJac(n_both, r);

	// compute sparsity pattern corresponding to Hessian w.r.t (theta, u)
	bool transpose = true;
	sparsity_pattern s(1), pattern;
	s[0] = true;
	pattern = a1_adfun_.RevSparseHes(n_both, s, transpose);

	// Determine row and column indices in lower triangle of Hessian w.r.t u
	// u indices start at n_fixed_
	for(size_t i = n_fixed_; i < n_both; i++)
	{	for(size_t j = n_fixed_; j < n_both; j++)
		{	// only compute lower triangular of Hessian w.r.t u only
			if( pattern[i * n_both + j] && i >= j )
			{	// only compute lower triangular of Hessian w.r.t u only
				if( i >= j )
				{	row_.push_back(i);
					col_.push_back(j);
				}
			}
		}
	}

	// create a weighting vector
	a1d_vector a1_w(1);
	a1_w[0] = 1.0;

	// place where results go (not used here)
	a1d_vector a1_val_out( row_.size() );

	// compute the work vector
	CppAD::sparse_hessian_work work;
	a1_adfun_.SparseHessian(
		a1_theta_u,
		a1_w,
		pattern,
		row_,
		col_,
		a1_val_out,
		work_
	);
}
//
// operator()
void newton_step_algo::operator()(
	const a1d_vector& a1_theta_u_v    ,
	a1d_vector&       a1_logdet_step  )
{	assert( a1_theta_u_v.size() == n_fixed_ + 2 * n_random_ );
	assert( a1_logdet_step.size() == 1 + n_random_ );

	// extract theta_u subvector
	a1d_vector a1_theta_u(n_fixed_ + n_random_);
	for(size_t j = 0; j < n_fixed_ + n_random_; j++)
		a1_theta_u[j] = a1_theta_u_v[j];

	// sparsity pattern not needed once we have work_
	CppAD::vectorBool not_used;

	// compute the sparse Hessian
	a1d_vector a1_w(1), a1_val_out( row_.size() );
	a1_w[0] = 1.0;
	a1_adfun_.SparseHessian(
		a1_theta_u,
		a1_w,
		not_used,
		row_,
		col_,
		a1_val_out,
		work_
	);

	// declare eigen matrix types
	using Eigen::Dynamic;
	typedef Eigen::Matrix<a1_double, Dynamic, Dynamic> dense_matrix;
	typedef Eigen::SparseMatrix<a1_double>             sparse_matrix;

	// create a lower triangular eigen sparse matrix representation of Hessian
	sparse_matrix hessian(n_random_, n_random_);
	size_t K = row_.size();
	for(size_t k = 0; k < K; k++)
	{	assert( n_fixed_ <= col_[k]  );
		assert( col_[k]  <= row_[k] );
		size_t i = row_[k] - n_fixed_;
		size_t j = col_[k] - n_fixed_;
		hessian.insert(i, j) = a1_val_out[k];
	}

	// compute an LDL^T Cholesky factorization of f_{u,u}(theta, u)
	Eigen::SimplicialLDLT<sparse_matrix, Eigen::Lower> chol;
	chol.analyzePattern(hessian);
	chol.factorize(hessian);

	// solve the equation f_{uu}{(2)}( theta, u ) * step = v
	dense_matrix rhs(n_random_, 1), step(n_random_, 1);
	for(size_t j = 0; j < n_random_; j++)
		rhs(j) = a1_theta_u_v[n_fixed_ + n_random_ + j];
	step = chol.solve(rhs);

	// compute the logdet( f_{u,u}(theta, u )
	dense_matrix diag = chol.vectorD();
	a1_double logdet = a1_double(0.0);
	for(size_t j = 0; j < n_random_; j++)
		logdet += log( diag(j) );

	// pack resutls into return vector
	a1_logdet_step[0] = logdet;
	for(size_t j = 0; j < n_random_; j++)
		a1_logdet_step[1 + j] = step(j);

	return;
}
// -------------------------------------------------------------------------
// newton_step ctor
newton_step::newton_step(void)
: atom_fun_(CPPAD_MIXED_NULL_PTR)
{ }
// newton_step destructor
newton_step::~newton_step(void)
{	if( atom_fun_ != CPPAD_MIXED_NULL_PTR )
		delete atom_fun_;
}
// initialize
void newton_step::initialize(
	CppAD::ADFun<a1_double>&          a1_adfun   ,
	const CppAD::vector<double>&      theta      ,
	const CppAD::vector<double>&      u          )
{	assert( atom_fun_ == CPPAD_MIXED_NULL_PTR );
	//
	size_t n_fixed  = theta.size();
	size_t n_random = u.size();
	// algo
	newton_step_algo algo(a1_adfun, theta, u);
	// atom_fun_
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
	atom_fun_ = new CppAD::checkpoint<double>(
		name, algo, a1_theta_u_v, a1_logdet_step, sparsity
	);
	assert( atom_fun_ != CPPAD_MIXED_NULL_PTR );
}
// size_var
size_t newton_step::size_var(void)
{	if( atom_fun_ ==  CPPAD_MIXED_NULL_PTR )
		return 0;
	return atom_fun_->size_var();
}
// eval
void newton_step::eval(
	const a1d_vector& a1_theta_u_v  ,
	a1d_vector&       a1_logdet_step )
{	assert( atom_fun_ != CPPAD_MIXED_NULL_PTR );
	(*atom_fun_)(a1_theta_u_v, a1_logdet_step);
}

} } // END_CPPAD_MIXED_NAMESPACE
