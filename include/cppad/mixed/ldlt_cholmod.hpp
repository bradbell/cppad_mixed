// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_LDLT_CHOLMOD_HPP
# define CPPAD_MIXED_LDLT_CHOLMOD_HPP

/*
$begin ldlt_cholmod$$
$spell
	suitesparse
	suitesparse
	Cholmod
	CppAD
	cholesky
	chol
	hes
$$

$section A Cholmod Cholesky Factor Class$$

$head See Also$$
$cref ldlt_eigen$$

$head Private$$
This class is an implementation detail and not part of the
CppAD Mixed user API.

$head Purpose$$
This class has utilities that work with a $code cholmod$$ LDLT factor,
$code cholmod$$ is part of the
$cref/suitesparse/install_unix/System Requirements/suitesparse/$$ package.

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


$head Example$$
The file $cref ldlt_cholmod.cpp$$ contains an example and test
using the operations in this class.

$head Preprocessor Symbols$$
The following extra $codei%CHOLMOD_%*%$$ symbols are defined
$srcthisfile%4%// BEGIN SYMBOLS%// END SYMBOLS%1%$$


$childtable%src/cholmod/constructor.cpp
	%src/cholmod/init.cpp
	%src/cholmod/pattern.cpp
	%src/cholmod/update.cpp
	%src/cholmod/logdet.cpp
	%src/cholmod/solve_H.cpp
	%src/cholmod/sim_cov.cpp
	%src/cholmod/inv.cpp
	%example/private/ldlt_cholmod.cpp
	%example/private/cholmod_factor.cpp
	%example/private/cholmod_solve.cpp
	%example/private/cholmod_solve2_a.cpp
	%example/private/cholmod_solve2_sim.cpp
%$$

$end
------------------------------------------------------------------------------
*/

# include <cppad/mixed/include_cholmod.hpp>
# include <cppad/mixed/typedef.hpp>

// BEGIN SYMBOLS
# define CHOLMOD_TRUE                  1
# define CHOLMOD_FALSE                 0
# define CHOLMOD_STYPE_NOT_SYMMETRIC   0
# define CHOLMOD_STYPE_LOWER_TRIANGLE -1
// END SYMBOLS


namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

class ldlt_cholmod {
private:
	//
	const size_t    nrow_;           // number of rows (and columns) in H
	bool            init_done_;     // has init been called
	bool            update_called_; // has update been called
	sparse_rc       H_rc_;          // sparsity pattern for H
	//
	CppAD::vector<size_t> key_;    // temporary sorting keys
	CppAD::vector<size_t> index_;  // temporary sorting indices
	//
	// mapping from sparse_mat_info to cholmod_sparse order
	CppAD::vector<size_t> H_rc2cholmod_order_;
	//
	cholmod_common    common_;
	cholmod_sparse*   sym_matrix_; // The symmetric matrix we are factoring
	cholmod_factor*   factor_;     // lower triangular LDL' factor
	//
	// information used by inv and computed during first call to inv
	CppAD::vector<size_t> sparseinv_order_;
	CppAD::vector<int>    sparseinv_p_, sparseinv_i_;
	//
	//
	// work space for solving equations
	cholmod_dense*    rhs_;        // right hand side of equation
	cholmod_sparse*   rhs_set_;    // sparsity pattern for rhs
	cholmod_dense*    sol_;        // solution of equaiton
	cholmod_sparse*   sol_set_;    // sparsity pattern for solutions
	cholmod_dense*    work_one_;   // first work space for solving equations
	cholmod_dense*    work_two_;   // second work space for solving equations
public:
	// ----------------------------------------------------------------------
	// non-const functions
	//
	// constructor
	ldlt_cholmod(size_t nrow);
	// destructor
	~ldlt_cholmod(void);
	// initialize
	void init(const sparse_rc& hes_rc);
	// factorize
	bool update(const d_sparse_rcv& hes_rcv);
	// ----------------------------------------------------------------------
	// const functions
	//
	// pattern
	const sparse_rc& pattern(void) const;
	//
	// log determinant
	double logdet(size_t& negative) const;
	//
	// compute a subset of the inverse
	void inv(
		const CppAD::vector<size_t>& row      ,
		const CppAD::vector<size_t>& col      ,
		CppAD::vector<double>&       val
	);
	// solve linear equations
	void solve_H(
		const CppAD::vector<size_t>& row      ,
		const CppAD::vector<double>& val_in   ,
		CppAD::vector<double>&       val_out
	);
	// simualte covariance
	bool sim_cov(
		const CppAD::vector<double>& w ,
		CppAD::vector<double>&       v
	);
};


} } // END_CPPAD_MIXED_NAMESPACE


# endif
