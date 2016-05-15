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

$head See Also$$
$cref ldlt_cholmod$$

$head Private$$
This class is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Factorization$$
The factorization is
$latex \[
	L D L^\R{T} = P H P^{T}
\] $$
where

$subhead H$$
is the matrix corresponding the current
$cref/update/ldlt_cholmod_update/$$.

$subhead L$$
is a lower triangular matrix with ones on the diagonal,

$subhead D$$
is a diagonal matrix.

$subhead P$$
is a permutation matrix.

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

$head Example$$
The file $cref ldlt_eigen_xam.cpp$$ contains an example and test
using the operations in this class.


$childtable%src/eigen/ldlt_eigen.cpp
	%example/private/ldlt_eigen_xam.cpp
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
	const size_t    n_row_;
	eigen_ldlt*     ptr_;
	//
public:
	// constructor
	ldlt_eigen(size_t n_row);
	// destructor
	~ldlt_eigen(void);
	// init
	void init(const CppAD::mixed::sparse_mat_info& hes_info);
	// update
	bool update(const CppAD::mixed::sparse_mat_info& hes_info);
	// logdet
	double logdet(size_t& negative) const;
	// solve
	void solve_H(
		const CppAD::vector<size_t>& row     ,
		const CppAD::vector<double>& val_in  ,
		CppAD::vector<double>&       val_out
	);
	// sim_cov
	bool sim_cov(
	const CppAD::vector<double>& w  ,
	CppAD::vector<double>&       v
	);
};


} } // END_CPPAD_MIXED_NAMESPACE


# endif
