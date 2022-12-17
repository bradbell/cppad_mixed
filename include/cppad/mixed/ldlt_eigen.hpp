/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
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

$nospell
$bold This is old cppad_mixed documentation:$$ Here is a link to its
$href%http://bradbell.github.io/cppad_mixed%current documentation%$$.
$$
$section An Eigen LDLT Factor Class$$

$head See Also$$
$cref ldlt_cholmod$$

$head Private$$
This class is an implementation detail and not part of the
CppAD Mixed user API.

$head Factorization$$
The factorization is
$latex \[
	L D L^\R{T} = P H P^{T}
\] $$

$subhead H$$
is the matrix corresponding the current
$cref/update/ldlt_cholmod_update/$$.

$subhead L$$
is a lower triangular matrix with ones on the diagonal,

$subhead D$$
is a diagonal matrix.

$subhead P$$
is a permutation matrix.

$head Double$$
This is the type of the elements in the matrices and is either
 double29628 or  a1_double29628.

$head Example$$
The file $cref ldlt_eigen.cpp$$ contains an example and test
using the operations in this class.


$childtable%src/eigen/ldlt_eigen.cpp
	%example/private/ldlt_eigen.cpp
%$$

$end
------------------------------------------------------------------------------
*/

# include <Eigen/Sparse>
# include <cppad/cppad.hpp>
# include <cppad/mixed/typedef.hpp>


namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

template <typename Double>
class ldlt_eigen {
private:
    typedef CppAD::vector<Double>                             v_vector;
	typedef Eigen::SparseMatrix<Double, Eigen::ColMajor>      eigen_sparse;
	typedef Eigen::Matrix<Double, Eigen::Dynamic, 1>          eigen_vector;
	typedef Eigen::PermutationMatrix<Eigen::Dynamic>          eigen_perm;
	typedef Eigen::SimplicialLDLT<eigen_sparse, Eigen::Lower> eigen_ldlt;
	//
	const size_t    n_row_;         // number of rows (and columns) in H
	bool            init_done_;     // has init been called
	bool            update_called_; // has update been called
	sparse_rc       H_rc_;          // sparsity pattern for H
	//
	eigen_ldlt*     ptr_;           // eigens ldlt factorization
	//
public:
	// ----------------------------------------------------------------------
	// non-const functions
	//
	// constructor
	ldlt_eigen(size_t n_row);
	// destructor
	~ldlt_eigen(void);
	// init
	void init(const sparse_rc& hes_rc);
	// update
	bool update(const CppAD::sparse_rcv<s_vector, v_vector>& hes_rcv);
	// ----------------------------------------------------------------------
	// const functions
	//
	// pattern
	const sparse_rc& pattern(void) const;
	//
	// split
	void split(
		eigen_sparse& L ,
		eigen_vector& D ,
		eigen_perm&   P
	) const;
	//
	// logdet
	Double logdet(size_t& negative) const;
	// solve
	void solve_H(
		const CppAD::vector<size_t>& row     ,
		const CppAD::vector<Double>& val_in  ,
		CppAD::vector<Double>&       val_out
	) const;
	//
	// sim_cov
	bool sim_cov(
	const CppAD::vector<Double>& w  ,
	CppAD::vector<Double>&       v
	) const;
	//
	// compute a subset of the inverse
	void inv(
		const CppAD::vector<size_t>& row      ,
		const CppAD::vector<size_t>& col      ,
		CppAD::vector<Double>&       val
	) const;
	// ----------------------------------------------------------------------
	// static functions
	//
	// solve_LDLT
	static eigen_vector solve_LDLT(
		const eigen_sparse& L  ,
		const eigen_vector& D  ,
		const eigen_perm&   P  ,
		const eigen_vector& b
	);
};


} } // END_CPPAD_MIXED_NAMESPACE


# endif
