// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_SPARSE_UP_TRI_SOL_HPP
# define CPPAD_MIXED_SPARSE_UP_TRI_SOL_HPP

# include <Eigen/SparseCore>

namespace CppAD { namespace mixed {
	Eigen::SparseMatrix<double, Eigen::ColMajor> sparse_up_tri_sol(
		const Eigen::SparseMatrix<double, Eigen::RowMajor>&  left  ,
		const Eigen::SparseMatrix<double, Eigen::ColMajor>&  right
	);
} }

# endif
