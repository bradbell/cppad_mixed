// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_LDLT_CHOLMOD_HPP
# define CPPAD_MIXED_LDLT_CHOLMOD_HPP

/*
$begin ldlt_cholmod$$
$spell
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
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
This class has utilities that work with a $code cholmod$$ LDLT factor,
$code cholmod$$ is part of the
$cref/SuiteSparse/install_unix/Special Requirements/SuiteSparse/$$ package.

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
The file $cref ldlt_cholmod_xam.cpp$$ contains an example and test
using the operations in this class.

$head Preprocessor Symbols$$
The following extra $codei%CHOLMOD_%*%$$ symbols are defined
$srcfile%include/cppad/mixed/ldlt_cholmod.hpp
	%4%// BEGIN SYMBOLS%// END SYMBOLS%1%$$


$childtable%src/cholmod/constructor.cpp
	%src/cholmod/init.cpp
	%src/cholmod/update.cpp
	%src/cholmod/logdet.cpp
	%src/cholmod/solve_H.cpp
	%src/cholmod/sim_cov.cpp
	%example/private/ldlt_cholmod_xam.cpp
	%example/private/cholmod_solve_xam.cpp
	%example/private/cholmod_solve2_a.cpp
	%example/private/cholmod_solve2_sim.cpp
%$$

$end
------------------------------------------------------------------------------
*/

# include <cholmod.h>
# include <cppad/mixed/sparse_mat_info.hpp>

// BEGIN SYMBOLS
# define CHOLMOD_TRUE                  1
# define CHOLMOD_FALSE                 0
# define CHOLMOD_STYPE_NOT_SYMMETRIC   0
# define CHOLMOD_STYPE_LOWER_TRIANGLE -1
// END SYMBOLS


namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

class ldlt_cholmod {
private:
	const size_t          nrow_;   // number of rows in sym_matrix_
	CppAD::vector<size_t> key_;    // temporary sorting keys
	CppAD::vector<size_t> index_;  // temporary sorting indices
	// mapping from sparse_mat_info to cholmod_sparse order
	CppAD::vector<size_t> info2cholmod_order_;
	//
	cholmod_common    common_;
	cholmod_sparse*   sym_matrix_; // The symmetric matrix we are factoring
	cholmod_factor*   factor_;     // lower triangular LDL' factor
	cholmod_dense*    rhs_;        // right hand side of equation
	cholmod_sparse*   rhs_set_;    // sparsity pattern for rhs
	cholmod_dense*    sol_;        // solution of equaiton
	cholmod_sparse*   sol_set_;    // sparsity pattern for solutions
	cholmod_dense*    work_one_;   // first work space for solving equations
	cholmod_dense*    work_two_;   // second work space for solving equations
public:
	// constructor
	ldlt_cholmod(size_t nrow);
	// destructor
	~ldlt_cholmod(void);
	// initialize
	void init(const sparse_mat_info& hes_info);
	// factorize
	bool update(const sparse_mat_info& hes_info);
	// log determinant
	double logdet(size_t& negative) const;
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
