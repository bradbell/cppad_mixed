// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/cppad_mixed.hpp>

/*
$begin init_fix_like$$
$spell
	CppAD
	init
	cppad
	jac
	vec
	const
	Cpp
	Jacobian
	var
	hes
$$

$section Initialize Fixed Likelihood$$

$head Syntax$$
$icode%mixed_object%.init_fix_like(%fixed_vec%)%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head fixed_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fixed_vec%
%$$
It specifies the value of the
$cref/fixed effects/cppad_mixed/Fixed Effects, theta/$$
vector $latex \theta$$ at which the initialization is done.

$head init_fix_like_done_$$
When this function is called, this member variable must be false.
Upon return it is true.

$head fix_like_fun_$$
On input, the member variable
$codei%
	CppAD::ADFun<double> fix_like_fun_
%$$
must be empty; i.e., $code fix_like_fun_.size_var() == 0$$.
If the return value for
$cref/fix_like/fix_like/$$ is empty,
$code fix_like_fun_$$ is not modified.
Otherwise,
upon return it contains the corresponding recording for the
$cref/fix_like/fix_like/$$.
The function result is the
$cref/negative log-density vector/cppad_mixed/Negative Log-Density Vector/$$
corresponding to the function
$cref/g(theta)/theory/Fixed Likelihood, g(theta)/$$.

$head fix_like_jac_$$
The input value of
$codei%
	sparse_jac_info fix_like_jac_
%$$
does not matter.
If $icode quasi_fixed$$ is false,
upon return $code fix_like_jac_$$ contains
$cref sparse_jac_info$$ for the
Jacobian corresponding to
$latex g^{(1)}) ( \theta )$$ see
$cref/g(theta)/theory/Fixed Likelihood, g(theta)/$$.

$subhead fix_like_fun_$$
This ADFun object can be used for the
$cref/sparse Jacobian call/sparse_jac_info/Sparse Jacobian Call/f/$$.

$head fix_like_hes_$$
The input value of
$codei%
	sparse_hes_info fix_like_hes_
%$$
does not matter.
If $icode quasi_fixed$$ is false,
upon return $code fix_like_hes_$$ contains
$cref sparse_hes_info$$ for the
lower triangle of a Hessian corresponding to
$latex g^{(2)}) ( \theta )$$ see
$cref/g(theta)/theory/Fixed Likelihood, g(theta)/$$.

$subhead fix_like_fun_$$
This ADFun object can be used for the
$cref/sparse Hessian call/sparse_hes_info/Sparse Hessian Call/f/$$.

$end
*/
# define DEBUG_FIX_LIKE_FUN 0


void cppad_mixed::init_fix_like(const d_vector& fixed_vec  )
{	assert( fixed_vec.size() == n_fixed_ );
	assert( ! init_fix_like_done_ );

	// ------------------------------------------------------------------------
	// fix_like_fun_
	// ------------------------------------------------------------------------
	// convert to an a1d_vector
	a1d_vector a1_theta(n_fixed_);
	for(size_t j = 0; j < n_fixed_; j++)
		a1_theta[j] = fixed_vec[j];

	// start recording a1_double operations
# if DEBUG_FIX_LIKE_FUN
	{	assert( n_fixed_ == 98 );
		assert( n_fixed_ * sizeof(double) == 784 );
		const char* file_name = "/tmp/fileJBRDSl";
		std::ifstream infile(file_name, std::ifstream::binary);
		CppAD::vector<double> double_buffer(n_fixed_);
		char   *char_buffer   = reinterpret_cast<char *>(double_buffer.data());
		infile.read(char_buffer, n_fixed_ * sizeof(double) );
		for(size_t j = 0; j < n_fixed_; j++)
			a1_theta[j] = double_buffer[j];
	}
# endif
	Independent(a1_theta);

	// compute fix_like
	a1d_vector a1_vec = fix_like(a1_theta);
	if( a1_vec.size() == 0 )
	{	CppAD::AD<double>::abort_recording();
		init_fix_like_done_ = true;
		assert( fix_like_fun_.size_var() == 0 );
		return;
	}

	// save the recording
	fix_like_fun_.Dependent(a1_theta, a1_vec);

	// optimize the recording
# ifdef NDEBUG
	fix_like_fun_.optimize();
# endif

	// ------------------------------------------------------------------------
	// fix_like_jac_.row, fix_like_jac_.col, fix_like_jac_.work
	// ------------------------------------------------------------------------
	// compute the sparsity pattern for the Jacobian
	typedef CppAD::vector< std::set<size_t> > sparsity_pattern;
	sparsity_pattern r(n_fixed_);
	for(size_t j = 0; j < n_fixed_; j++)
		r[j].insert(j);
	sparsity_pattern pattern = fix_like_fun_.ForSparseJac(n_fixed_, r);

	// convert sparsity to row and column index form
	fix_like_jac_.row.clear();
	fix_like_jac_.col.clear();
	std::set<size_t>::iterator itr;
	for(size_t i = 0; i < fix_like_fun_.Range(); i++)
	{	for(itr = pattern[i].begin(); itr != pattern[i].end(); itr++)
		{	size_t j = *itr;
			fix_like_jac_.row.push_back(i);
			fix_like_jac_.col.push_back(j);
		}
	}

	// direction for this sparse Jacobian
	fix_like_jac_.direction = sparse_jac_info::Forward;

	// compute the work vector for reuse during Jacobian sparsity calculations
	d_vector jac( fix_like_jac_.row.size() );
	fix_like_fun_.SparseJacobianForward(
		fixed_vec       ,
		pattern         ,
		fix_like_jac_.row  ,
		fix_like_jac_.col  ,
		jac             ,
		fix_like_jac_.work
	);
	if( quasi_fixed_ )
	{	init_fix_like_done_ = true;
		return;
	}
	// ------------------------------------------------------------------------
	// fix_like_hes_.row, fix_like_hes_.col, fix_like_hes_.work
	// ------------------------------------------------------------------------
	// no need to recalculate forward sparsity pattern.
	//
	// sparsity pattern for the Hessian
	sparsity_pattern s(1);
	for(size_t i = 0; i < fix_like_fun_.Range(); i++ )
		s[0].insert(i);
	pattern.clear();
	pattern = fix_like_fun_.RevSparseHes(n_fixed_, s);

	// determine row and column indices in lower triangle of Hessian
	fix_like_hes_.row.clear();
	fix_like_hes_.col.clear();
	for(size_t i = 0; i < n_fixed_; i++)
	{	for(itr = pattern[i].begin(); itr != pattern[i].end(); itr++)
		{	size_t j = *itr;
			// only compute lower triangular part
			if( i >= j )
			{	fix_like_hes_.row.push_back(i);
				fix_like_hes_.col.push_back(j);
			}
		}
	}
	size_t K = fix_like_hes_.row.size();

	// compute the work vector for reuse during Hessian sparsity calculations
	d_vector weight( fix_like_fun_.Range() ), hes(K);
	fix_like_fun_.SparseHessian(
		fixed_vec       ,
		weight          ,
		pattern         ,
		fix_like_hes_.row  ,
		fix_like_hes_.col  ,
		hes             ,
		fix_like_hes_.work
	);

	init_fix_like_done_ = true;
	return;
}


