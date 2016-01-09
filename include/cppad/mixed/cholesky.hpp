// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_CHOLESKY_HPP
# define CPPAD_MIXED_CHOLESKY_HPP

/*
$begin cholesky$$
$spell
	CppAD
	cholesky
	chol
	hes
$$

$section Cholesky Factor Class$$

$head Syntax$$
$icode%CppAD::mixed::cholesky %chol_hes_ran%%;
%$$

$head chol_hes_ran$$
This constructor creates a Cholesky factor.
It is called $icode chol_hes_ran$$ because it is only intended for
a cholesky factor of the Hessian with respect to the random effects; i.e.,
$latex f_{u,u} ( \theta , u )$$.

$childtable%src/eigen/cholesky.cpp
%$$

$end
------------------------------------------------------------------------------
*/

# include <Eigen/Sparse>
# include <cppad/cppad.hpp>


namespace Eigen {
	// template classes
	template <typename _Scalar, int _Options, typename _Index>
	class SparseMatrix;

	template <typename Index>
	class AMDOrdering;

	template <typename _MatrixType, int _UpLo, typename _Ordering>
	class SimplicialLDLT;
}

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

class cholesky {
public:
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor>      eigen_sparse;
	typedef Eigen::SimplicialLDLT<eigen_sparse, Eigen::Lower> eigen_cholesky;
private:
	eigen_cholesky* ptr_;
	//
public:
	// constructor
	cholesky(void);
	// destructor
	~cholesky(void);
	// ptr
	const eigen_cholesky* ptr(void) const
	{	return ptr_; }
	// analyze
	void analyze(
		size_t                       n_fixed  ,
		size_t                       n_random ,
		const CppAD::vector<size_t>& row      ,
		const CppAD::vector<size_t>& col
	);
	// factorize
	void factorize(
		size_t                       n_fixed      ,
		size_t                       n_random     ,
		const CppAD::vector<size_t>& row          ,
		const CppAD::vector<size_t>& col          ,
		const CppAD::vector<double>& both         ,
		CppAD::ADFun<double>&        hes_ran_fun
	);
	// logdet
	double logdet(size_t n_random) const;
};


} } // END_CPPAD_MIXED_NAMESPACE


# endif
