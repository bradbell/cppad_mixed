// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_SPARSE_JAC_INFO_HPP
# define CPPAD_MIXED_SPARSE_JAC_INFO_HPP
/*
$begin sparse_jac_info$$
$spell
	CppAD
	const
	std
	bool
	work work
	jac
	CppAD
	Jacobian
%$$

$section Sparse Jacobian Information$$

$head Syntax$$
$codei%sparse_jac_info %jac_info%$$

$head Purpose$$
This structure holds information about a specific Jacobian.

$head row$$
The field $icode%jac_info%.row%$$ has prototype
$codei%
	CppAD::vector<size_t> %jac_info%.row
%$$
It has size zero when it is constructed.
After initialization it should contain the row indices
corresponding to possibly non-zero elements of the Jacobian.

$head col$$
The field $icode%jac_info%.col%$$ has prototype
$codei%
	CppAD::vector<size_t> %jac_info%.col
%$$
It has size zero when it is constructed.
After initialization it should contain the column indices
corresponding to possibly non-zero elements of the Jacobian
and has the same size as $icode%jac_info%.row%$$.

$head direction$$
The field $icode%jac_info%.direction%$$ has prototype
$codei%
	Direction %direction%
%$$
and its value is either
$code sparse_jac_info::Forward$$ or
$code sparse_jac_info::Reverse$$.
It determines the value of $icode Direction$$
($code Forward$$ or $code Reverse$$)
in the sparse Jacobian call.

$head work$$
The field $icode%jac_info%.work%$$ has prototype
$codei%
	CppAD::sparse_jacobian_work work
%$$
It has no information when it is constructed.
After initialization it should contain the CppAD cached information.

$head Sparse Jacobian Call$$
$codei%
	%f%.SparseJacobian%Direction%(
		%x%,
		%not_used%,
		%jac_info%.row,
		%jac_info%.col,
		%jac%,
		%jac_info%.work
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
It is the location where the Jacobian is being evaluated.

$subhead w$$
The argument $icode w%$$ has prototype
$codei%
	const CppAD::vector<%Base%>& %w%
%$$
and its size is $icode%f%.Range()%$$.
It is the weighting for the components of the Jacobian that is
being computed.

$subhead not_used$$
This argument has one of the following prototypes
$codei%
	CppAD::vector < std::set<size_t> > %not_used%
	CppAD::vector<bool>  %not_used%
	CppAD::vectorBool    %not_used%
%$$
It does not matter which and it is not used.

$subhead jac$$
This argument has prototype
$icode%
	CppAD::vector<%Base%>& %jac%
%$$
Its input size must be the same as $icode%jac_info%.row.size()%$$
and the value of its elements does not matter.
Upon return $icode%jac%[%k%]%$$ is the Jacobian at row index
$icode%jac_info%.row[%k%]%$$ and column index
$icode%jac_info%.col[%k%]%$$.


$end
*/

	struct sparse_jac_info {
		enum Direction { Forward , Reverse } direction;
		CppAD::vector<size_t>                row;
		CppAD::vector<size_t>                col;
		CppAD::sparse_jacobian_work          work;
	};

# endif
