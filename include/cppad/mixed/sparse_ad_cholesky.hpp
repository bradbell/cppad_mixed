// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_SPARSE_AD_CHOLESKY_HPP
# define CPPAD_MIXED_SPARSE_AD_CHOLESKY_HPP
/*
$begin sparse_ad_cholesky$$
$spell
	CppAD
	Cholesky
$$

$section Sparse Cholesky Factorization as an Atomic CppAD Operation$$

$head Purpose$$
Given a symmetric positive definite matrix $latex A \in \B{R}^{n \times n}$$,
this routine computes a permutation matrix $latex P \in \B{R}^{n \times n}$$
and a lower triangular matrix $latex L \in \B{R}^{n \times n}$$
such that
$latex \[
	P A P^\R{T} = L L^\R{T}
\] $$
The permutation matrix is chosen to be fill reducing; i.e.,
to make $latex L$$ sparse.

$head Notation$$

$subhead A$$
We use $latex A$$ to refer to the symmetric positive definite
matrix.

$subhead P$$
We use $latex P$$ to refer to the permutation matrix.

$subhead L$$
We use $latex L$$ to refer to the lower triangular matrix.

$head Private$$
This class is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Public Member Functions$$
$srcfile%include/cppad/mixed/sparse_ad_cholesky.hpp
	%4%// BEGIN PUBLIC MEMBER FUNCTIONS%// END PUBLIC MEMBER FUNCTIONS%1%$$

$head Type Declarations$$
$srcfile%include/cppad/mixed/sparse_ad_cholesky.hpp
	%4%// BEGIN TYPE DECLARATIONS%// END TYPE DECLARATIONS%1%$$

$head Member Variables$$
$srcfile%include/cppad/mixed/sparse_ad_cholesky.hpp
	%4%// BEGIN MEMBER VARIABLES%// END MEMBER VARIABLES%1%$$

$childtable%src/eigen/sparse_ad_cholesky.cpp
	%example/private/sparse_ad_cholesky_ad.cpp
	%example/private/sparse_ad_cholesky_p.cpp
	%example/private/sparse_ad_cholesky_eq.cpp
	%example/private/sparse_ad_cholesky_v.cpp
%$$

$end
*/
# include <cppad/cppad.hpp>
# include <Eigen/SparseCore>
# include <Eigen/SparseCholesky>
# include <cppad/mixed/sparse_mat_info.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// -----------------------------------------------------------------
// define AD, Dynamic, sparse_mat2low, some vector and matrix types
// -----------------------------------------------------------------
//

class sparse_ad_cholesky : public CppAD::atomic_base<double> {
// -----------------------------------------------------------------
// BEGIN TYPE DECLARATIONS
private:
	typedef Eigen::
	SparseMatrix<double, Eigen::ColMajor>               sparse_d_matrix;
	typedef Eigen::
	SparseMatrix< CppAD::AD<double>, Eigen::ColMajor >  sparse_ad_matrix;
// END TYPE DECLARATIONS
// -----------------------------------------------------------------
// BEGIN MEMBER VARIABLES
private:
	// number of columns and rows in the square matrices Alow and L
	const size_t nc_;
	//
	// OK flag is initialized as true by constructor and set to false if an
	// error occurs.
	bool ok_;
	//
	// Sparsity pattern for Alow and the temporary vector Alow_pattern_.val
	// (set by constructor).
	CppAD::mixed::sparse_mat_info Alow_pattern_;

	// Object used for Cholesky factorization
	// (analyzePattern is only called by constructor).
	Eigen::SimplicialLDLT<sparse_d_matrix> ldlt_obj_;

	// Value of the permutation matrix (set by constructor).
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P_;

	// Sparsity pattern for L and the temporary vector L_pattern_.val
	// (set by constructor).
	CppAD::mixed::sparse_mat_info L_pattern_;
	//
	// Indices that access lower traingle of B = P * A * P^T
	// in column major order (set by constructor).
	CppAD::vector<size_t> Alow_permuted_;
	//
	// Indices that access L_pattern_ in row major order
	// (set by constructor).
	CppAD::vector<size_t> L_row_major_;
// END MEMBER VARIABLES
// -----------------------------------------------------------------
// BEGIN PUBLIC MEMBER FUNCTIONS
public:
	//
	// OK status flag
	bool ok(void)
	{	return ok_; }
	//
	// constructor
	sparse_ad_cholesky(const sparse_d_matrix& Alow );
	//
	// Permutation matrix for this factorization
	const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>&
	permutation(void);
	//
	// AD version of atomic Cholesky factorization
	void eval(
		const sparse_ad_matrix& ad_Alow  ,
		sparse_ad_matrix&       ad_L
	);
// END PUBLIC MEMBER FUNCTIONS
// -----------------------------------------------------------------
// private virtual functions
private:
	// ------------------------------------------------------------------
	// CppAD forward mode for this operation (called by CppAD)
	virtual bool forward(
		// lowest order Taylor coefficient we are evaluating
		size_t                          p ,
		// highest order Taylor coefficient we are evaluating
		size_t                          q ,
		// which components of x are variables
		const CppAD::vector<bool>&      vx ,
		// which components of y are variables
		CppAD::vector<bool>&            vy ,
		// tx [ j * (q+1) + k ] is x_j^k
		const CppAD::vector<double>&    tx ,
		// ty [ i * (q+1) + k ] is y_i^k
		CppAD::vector<double>&          ty
	);
	// ------------------------------------------------------------------
	// reverse mode for this operation (called by CppAD)
	virtual bool reverse(
		// highest order Taylor coefficient that we are computing derivative of
		size_t                     q ,
		// forward mode Taylor coefficients for x variables
		const CppAD::vector<double>&     tx ,
		// forward mode Taylor coefficients for y variables
		const CppAD::vector<double>&     ty ,
		// upon return, derivative of G[ F[ {x_j^k} ] ] w.r.t {x_j^k}
		CppAD::vector<double>&           px ,
		// derivative of G[ {y_i^k} ] w.r.t. {y_i^k}
		const CppAD::vector<double>&     py
	);
}; // END_SPARSE_AD_CHOLESKY

} } // END_CPPAD_MIXED_NAMESPACE

# endif
