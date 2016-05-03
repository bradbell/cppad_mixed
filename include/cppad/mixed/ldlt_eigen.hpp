// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_LDLT_EIGEN_HPP
# define CPPAD_MIXED_LDLT_EIGEN_HPP

/*
$begin ldlt_eigen$$
$spell
	ldlt_eigen
	Simplicial
	CppAD
	cholesky
	chol
	hes
	eigen
	typedef
$$

$section An Eigen LDLT Factor Class$$

$head Private$$
This class is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head eigen_sparse$$
The type $code CppAD::mixed::ldlt_eigen::eigen_sparse$$ is defined by
$codei%
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor>  eigen_sparse;
%$$

$head eigen_ldlt_eigen$$
The type $code CppAD::mixed::ldlt_eigen::eigen_ldlt_eigen$$ is defined by
$codei%
	typedef Eigen::SimplicialLDLT<eigen_sparse, Eigen::Lower> eigen_ldlt;
%$$


$childtable%src/eigen/ldlt_eigen.cpp
%$$

$end
------------------------------------------------------------------------------
*/

# include <Eigen/Sparse>
# include <cppad/cppad.hpp>
# include <cppad/mixed/sparse_mat_info.hpp>


namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

class ldlt_eigen {
public:
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor>      eigen_sparse;
	typedef Eigen::SimplicialLDLT<eigen_sparse, Eigen::Lower> eigen_ldlt;
private:
	const size_t    n_random_;
	eigen_ldlt* ptr_;
	//
public:
	// constructor
	ldlt_eigen(size_t n_random);
	// destructor
	~ldlt_eigen(void);
	// init
	void init(const CppAD::mixed::sparse_mat_info& hes_info);
	// update
	bool update(const CppAD::mixed::sparse_mat_info& hes_info);
	// logdet
	double logdet(void) const;
	// solve
	void solve_H(
		const CppAD::vector<size_t>& row_in  ,
		const CppAD::vector<double>& val_in  ,
		CppAD::vector<size_t>&       row_out ,
		CppAD::vector<double>&       val_out
	);
};


} } // END_CPPAD_MIXED_NAMESPACE


# endif
