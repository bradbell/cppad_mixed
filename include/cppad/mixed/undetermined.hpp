// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
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
