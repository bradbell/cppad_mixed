// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_CHOL_RAN_HES_HPP
# define CPPAD_MIXED_CHOL_RAN_HES_HPP

# include <Eigen/Sparse>

namespace CppAD { namespace mixed {
   // This would be member of cppad_mixed class if it were not for all the
   // warnings Eigen generates. 2DO: review this because the Eigen warning
   // problem has been fixed using the SYSTEM include option.
   extern Eigen::SimplicialLDLT<
      Eigen::SparseMatrix<double> , Eigen::Lower
   > ldlt_ran_hes_;

   extern void analyze_ldlt_ran_hes(
      size_t                       n_fixed  ,
      size_t                       n_random ,
      const CppAD::vector<size_t>& row      ,
      const CppAD::vector<size_t>& col
   );
   extern void factorize_ldlt_ran_hes(
      size_t                       n_fixed  ,
      size_t                       n_random ,
      const CppAD::vector<size_t>& row      ,
      const CppAD::vector<size_t>& col      ,
      const CppAD::vector<double>& both     ,
      CppAD::ADFun<double>&        hessian
   );
   extern double logdet_ldlt_ran_hes(size_t n_random);
} }


# endif
