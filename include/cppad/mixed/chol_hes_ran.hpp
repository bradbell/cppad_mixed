// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_CHOL_HES_RAN_HPP
# define CPPAD_MIXED_CHOL_HES_RAN_HPP

# include <Eigen/Sparse>

namespace CppAD { namespace mixed {
	// This would be member of cppad_mixed class if it were not for all the
	// warnings Eigen generates.
	extern Eigen::SimplicialLDLT<
		Eigen::SparseMatrix<double> , Eigen::Lower
	> chol_hes_ran_;

	extern void analyze_chol_hes_ran(
		size_t                       n_fixed  ,
		size_t                       n_random ,
		const CppAD::vector<size_t>& row      ,
		const CppAD::vector<size_t>& col
	);
	extern void factorize_chol_hes_ran(
		size_t                       n_fixed  ,
		size_t                       n_random ,
		const CppAD::vector<size_t>& row      ,
		const CppAD::vector<size_t>& col      ,
		const CppAD::vector<double>& both     ,
		CppAD::ADFun<double>&        hessian
	);
	extern double logdet_chol_hes_ran(size_t n_random);
} }


# endif
