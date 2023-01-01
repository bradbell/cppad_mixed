/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin sparse_low2sym.cpp$$
$spell
	tri
	cppad
	sym
$$

$nospell
$bold This is cppad_mixed--20220519 documentation:$$ Here is a link to its
$href%https://cppad-mixed.readthedocs.io%current documentation%$$.
$$
$section sparse_low2sym: Example and Test$$

$head Private$$
This example is not part of the
$cref/cppad_mixed public API/base_class/$$.

$code
$srcthisfile%0%// BEGIN C++%// END C++%1%$$
$$
$end
*/
// BEGIN C++
# include <cmath>
# include <cstddef>
# include <limits>
# include <cppad/mixed/sparse_low2sym.hpp>

bool sparse_low2sym_xam(void)
{	bool ok = true;
	typedef Eigen::SparseMatrix<size_t> sparse_matrix;
	typedef Eigen::Index                Int;
	//
	size_t nr = 3;
	sparse_matrix lower;
	lower.resize(Int(nr), Int(nr));
	Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> check(nr, nr);
	size_t count = 0;
	for(size_t i = 0; i < nr; i++)
	{	for(size_t j = 0; j <= i; j++)
		{	lower.insert(Int(i), Int(j)) = count;
			check(i, j)                  = count;
			check(j, i)                  = count;
			++count;
		}
	}
	sparse_matrix symmetric = CppAD::mixed::sparse_low2sym(lower);
	//
	// check the result
	count = 0;
	for(size_t k = 0; k < nr; k++)
	{	for(sparse_matrix::InnerIterator itr(symmetric, Int(k)); itr; ++itr)
		{	Int i = itr.row();
			Int j = itr.col();
			ok   &= itr.value() == check(i, j);
			++count;
		}
	}
	ok &= count == nr * nr;
	return ok;
}
// END C++
