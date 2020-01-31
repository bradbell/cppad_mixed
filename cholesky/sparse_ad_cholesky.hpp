// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-20 University of Washington
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
	cppad
$$

$section Sparse Cholesky Factorization as an Atomic CppAD Operation$$

$head Syntax$$
$codei%CppAD::mixed::sparse_ad_cholesky %cholesky%$$

$head Not Used$$
The option for $code cppad_mixed$$ to use this class was removed on 2017-01-21.
Is implementation and examples have been stored in the
$code cholesky$$ directory for possible future research.

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
CppAD Mixed user API.

$head Public Member Functions$$
$srcthisfile%4%// BEGIN PUBLIC MEMBER FUNCTIONS%// END PUBLIC MEMBER FUNCTIONS%1%$$

$head Type Declarations$$
$srcthisfile%4%// BEGIN TYPE DECLARATIONS%// END TYPE DECLARATIONS%1%$$

$head Member Variables$$
$srcthisfile%4%// BEGIN MEMBER VARIABLES%// END MEMBER VARIABLES%1%$$

$childtable%cholesky/sparse_ad_cholesky.cpp
	%cholesky/example/sparse_ad_chol_eval.cpp
	%cholesky/example/sparse_ad_chol_perm.cpp
	%cholesky/example/sparse_ad_chol_eq.cpp
	%cholesky/example/sparse_ad_chol_var.cpp
	%cholesky/example/sparse_ad_chol_sp1.cpp
	%cholesky/example/sparse_ad_chol_sp2.cpp
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
    typedef CppAD::local::sparse::pack_setvec           pack_setvec;
	typedef Eigen::
	SparseMatrix<double, Eigen::ColMajor>               sparse_d_matrix;
	typedef Eigen::
	SparseMatrix< CppAD::AD<double>, Eigen::ColMajor >  sparse_ad_matrix;
// END TYPE DECLARATIONS
// -----------------------------------------------------------------
// BEGIN MEMBER VARIABLES
private:
	// number of columns and rows in the square matrices Alow and L
	// (set by initialize)
	size_t nc_;
	//
	// OK flag is initialized as true to false if an
	// error occurs.
	bool ok_;
	//
	// Sparsity pattern for Alow and the temporary vector Alow_pattern_.val
	// (set by initialize).
	CppAD::mixed::sparse_mat_info Alow_pattern_;

	// Object used for Cholesky factorization
	// (analyzePattern is only called by initialize).
	Eigen::SimplicialLDLT<sparse_d_matrix> ldlt_obj_;

	// Value of the permutation matrix (set by initialize).
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P_;

	// Sparsity pattern for L and the temporary vector L_pattern_.val
	// (set by initialize).
	CppAD::mixed::sparse_mat_info L_pattern_;
	//
	// Indices that access lower traingle of B = P * A * P^T
	// in column major order (set by initialize).
	CppAD::vector<size_t> Alow_permuted_;
	//
	// Indices that access L_pattern_ in row major order
	// (set by initialize).
	CppAD::vector<size_t> L_row_major_;
	//
	// Jacobian sparsity as a pack_setvec.  This is CppAD internal
	// representations for vectors of of bools.  Its use is not part of the
	// CppAD API (yet), but it is more efficient, so we use it here.
	pack_setvec jac_sparsity_pack_;
// END MEMBER VARIABLES
// -----------------------------------------------------------------
// BEGIN PUBLIC MEMBER FUNCTIONS
public:
	//
	// OK status flag
	bool ok(void)
	{	return ok_; }
	//
	// default constructor
	sparse_ad_cholesky(void) : CppAD::atomic_base<double>(
		"sparse_ad_cholesky",
		CppAD::atomic_base<double>::pack_sparsity_enum
	)
	{ }
	//
	// initialize
	void initialize(const sparse_ad_matrix& Alow );
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
// =======================================================================
private:
	// -------------------------------------------------------------------
	// private functions only used by this class
	// -------------------------------------------------------------------
	// compute the Jacobian sparsity patten for map Alow -> L
	template <class Sparsity> void set_jac_sparsity(
		Sparsity& jac_sparsity
	);
	//
	// compute the Hessian sparsity patten for map Alow -> L
	template <class Sparsity> void set_hes_sparsity(
		const CppAD::vector<bool>&    s            ,
		Sparsity&                     jac_sparsity ,
		Sparsity&                     hes_sparsity
	);
	//
	// ------------------------------------------------------------------
	// private virtual functions called and specified by CppAD
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
	// ------------------------------------------------------------------
	// vectorBool forward Jacobian sparsity pattern for S = f'(x) * R
	virtual bool for_sparse_jac(
		// number of columns in the matrix R
		size_t                        q         ,
		// sparsity pattern for R
		const CppAD::vectorBool&      r         ,
		// sparsity pattern for S
		CppAD::vectorBool&            s         ,
		// parameters in argument to atomic function
		const CppAD::vector<double>&  not_used
	);
	// ------------------------------------------------------------------
	// vectorBool reverse Jacobian sparsity pattern for S = R * f'(x)
	virtual bool rev_sparse_jac(
		// number of rows in the matrix R
		size_t                        q         ,
		// sparsity pattern for R^T
		const CppAD::vectorBool&      rt        ,
		// sparsity pattern for S^T
		CppAD::vectorBool&            st        ,
		// parameters in argument to atomic function
		const CppAD::vector<double>&  not_used
	);
	// ------------------------------------------------------------------
	// vectorBool reverse Hessian sparsity for V(x) = (g o f)^(2) (x) * R
	virtual bool rev_sparse_hes(
		// variable flag for x arguments
		const vector<bool>&                   vx       ,
		// sparsity pattern for scalar valued S(x) = g'[ f(x) ]
		const vector<bool>&                   s        ,
		// sparsity pattern for T(x) = (g o f)' (x) = S(x) * f'(x)
		vector<bool>&                         t        ,
		// number of columns in the matrix R
		size_t                                q        ,
		// sparsity pattern for the matrix R
		const CppAD::vectorBool&              r        ,
		// sparsity pattern for U(x) = g^(2)[ f(x) ] f'(x) R
		const CppAD::vectorBool&              u        ,
		// sparsity pattern for V(x) = (g o f)^(2) (x) * R
		CppAD::vectorBool&                    v        ,
		// parameters in argument to atomic function
		const vector<double>&                 not_used
	);
}; // END_SPARSE_AD_CHOLESKY

} } // END_CPPAD_MIXED_NAMESPACE

# endif
