// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_SPARSE_LOW_TRI_SOL_HPP
# define CPPAD_MIXED_SPARSE_LOW_TRI_SOL_HPP

# include <Eigen/SparseCore>

namespace CppAD { namespace mixed {
   Eigen::SparseMatrix<double, Eigen::ColMajor> sparse_low_tri_sol(
      const Eigen::SparseMatrix<double, Eigen::RowMajor>&  left  ,
      const Eigen::SparseMatrix<double, Eigen::ColMajor>&  right
   );
} }

# endif
