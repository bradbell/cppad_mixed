// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_UNDETERMINED_HPP
# define CPPAD_MIXED_UNDETERMINED_HPP

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

size_t undetermined(
	const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& A     ,
	const Eigen::Matrix<double, Eigen::Dynamic, 1>&              b     ,
	double                                                       delta ,
	Eigen::Matrix<size_t, Eigen::Dynamic, 1>&                    D     ,
	Eigen::Matrix<size_t, Eigen::Dynamic, 1>&                    I     ,
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&       C     ,
	Eigen::Matrix<double, Eigen::Dynamic, 1>&                    e
);

} } // END_CPPAD_MIXED_NAMESPACE

# endif
