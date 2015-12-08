// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_SPARSE_HES_INFO_HPP
# define CPPAD_MIXED_SPARSE_HES_INFO_HPP
/*
$begin sparse_hes_info$$
$spell
	CppAD
	const
	std
	bool
	work work
	hes
	CppAD
%$$

$section Sparse Hessian Information$$

$head Syntax$$
$codei%sparse_hes_info %hes_info%$$

$head Purpose$$
This structure holds information about a specific Hessian.

$head row$$
The field $icode%hes_info%.row%$$ has prototype
$codei%
	CppAD::vector<size_t> %hes_info%.row
%$$
It has size zero when it is constructed.
After initialization it should contain the row indices
corresponding to possibly non-zero elements of the Hessian.

$head col$$
The field $icode%hes_info%.col%$$ has prototype
$codei%
	CppAD::vector<size_t> %hes_info%.col
%$$
It has size zero when it is constructed.
After initialization it should contain the column indices
corresponding to possibly non-zero elements of the Hessian
and has the same size as $icode%hes_info%.row%$$.

$head work$$
The field $icode%hes_info%.work%$$ has prototype
$codei%
	CppAD::sparse_hessian_work work
%$$
It has no information when it is constructed.
After initialization it should contain the CppAD cached information.

$head Sparse Hessian Call$$
$codei%
	%f%.SparseHessian(
		%x%,
		%w%
		%not_used%,
		%hes_info%.row,
		%hes_info%.col,
		%hes%,
		%hes_info%.work
	)
%$$

$subhead f$$
The function $icode f$$ has prototype
$codei%
	CppAD::ADFun<%Base%> %f%
%$$

$subhead x$$
The argument $icode x%$$ has prototype
$codei%
	const CppAD::vector<%Base%>& %x%
%$$
and its size is $icode%f%.Domain()%$$.
It is the location where the Hessian is being evaluated.

$subhead w$$
The argument $icode w%$$ has prototype
$codei%
	const CppAD::vector<%Base%>& %w%
%$$
and its size is $icode%f%.Range()%$$.
It is the weighting for the components of the Hessian that is
being computed.

$subhead not_used$$
This argument has one of the following prototypes
$codei%
	CppAD::vector < std::set<size_t> > %not_used%
	CppAD::vector<bool>  %not_used%
	CppAD::vectorBool    %not_used%
%$$
It does not matter which and it is not used.

$subhead hes$$
This argument has prototype
$icode%
	CppAD::vector<%Base%>& %hes%
%$$
Its input size must be the same as $icode%hes_info%.row.size()%$$
and the value of its elements does not matter.
Upon return $icode%hes%[%k%]%$$ is the Hessian at row index
$icode%hes_info%.row[%k%]%$$ and column index
$icode%hes_info%.col[%k%]%$$.


$end
*/

	struct sparse_hes_info {
		CppAD::vector<size_t>      row;
		CppAD::vector<size_t>      col;
		CppAD::sparse_hessian_work work;
	};

# endif
