// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_CHOLMOD_HPP
# define CPPAD_MIXED_CHOLMOD_HPP

/*
$begin cholmod$$
$spell
	Cholmod
	CppAD
	cholesky
	chol
	hes
$$

$section Cholmod Cholesky Factor Class$$

$head Private$$
This class is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
This class has utilities that work with a $code cholmod$$ Cholesky factor
($code cholmod$$ is part of the
$cref/SuiteSparse/install_unix/Special Requirements/SuiteSparse/$$ package).
It is called $icode chol_ran_hes$$ because it is only intended for
the cholesky factor of the Hessian with respect to the random effects; i.e.,
$latex f_{u,u} ( \theta , u )$$.

$childtable%src/cholmod/constructor.cpp
	%src/cholmod/init.cpp
	%src/cholmod/update.cpp
%$$

$end
------------------------------------------------------------------------------
*/

# include <cholmod.h>
# include <cppad/mixed/sparse_mat_info.hpp>

# define CHOLMOD_TRUE                  1
# define CHOLMOD_FALSE                 0
# define CHOLMOD_STYPE_NOT_SYMMETRIC   0
# define CHOLMOD_STYPE_LOWER_TRIANGLE -1


namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

class cholmod {
private:
	const size_t      n_fixed_;    // number of fixed effects
	const size_t      n_random_;   // number of random effects
	cholmod_common    common_;
	cholmod_sparse*   pos_matrix_; // The positive matrix we are factoring
	cholmod_factor*   factor_;     // lower triangular LDL' factor
	cholmod_dense*    rhs_;        // right hand side of equation
	cholmod_sparse*   rhs_set_;    // sparsity pattern for rhs
	cholmod_dense*    sol_;        // solution of equaiton
	cholmod_sparse*   sol_set_;    // sparsity pattern for solutions
	cholmod_dense*    work_one_;   // first work space for solving equations
	cholmod_dense*    work_two_;   // second work space for solving equations
public:
	// constructor
	cholmod(size_t n_fixed, size_t n_random);
	// destructor
	~cholmod(void);
	// initialize
	void init(const sparse_mat_info& hes_info);
	// factorize
	void update(const sparse_mat_info& hes_info);
};


} } // END_CPPAD_MIXED_NAMESPACE


# endif
