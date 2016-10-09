// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/sparse_ad_cholesky.hpp>

/*
$begin sparse_ad_chol_sp.cpp$$
$spell
	Cholesky
$$

$section Sparse AD Cholesky Sparsity Calculation: Example and Test$$

$head Problem$$
We are given the function $latex A : \B{R}^3 \rightarrow \B{R}^{3 \times 3}$$
defined by
$latex \[
	A(x) = \left( \begin{array}{ccc}
		x_0 & 0    & x_2  \\
		0   & x_1  & 0   \\
		x_2 & 0    & x_3
	\end{array} \right)
\] $$

$head Permutation$$
The fill reducing permutation
$cref/P/sparse_ad_cholesky/Notation/P/$$
transposes indices zero and one; i.e.
$latex \[
	P= \left( \begin{array}{ccc}
		0   & 1    & 0    \\
		1   & 0    & 0   \\
		0   & 0    & 1
	\end{array} \right)
	\; , \;
	P A(x) P^\R{T} = \left( \begin{array}{ccc}
		x_1 & 0    & 0    \\
		0   & x_0  & x_2 \\
		0   & x_2  & x_3
	\end{array} \right)
\] $$

$head Cholesky Factor$$
The Cholesky factor
$cref/L/sparse_ad_cholesky/Notation/L/$$ is
$latex \[
	L(x) = \left( \begin{array}{ccc}
		\sqrt{x_1} & 0                   & 0    \\
		0          & \sqrt{x_0}          & 0   \\
		0          & x_2 / \sqrt{x_0}    & \sqrt{ x_3 - x_2^2 / x_0 }
	\end{array} \right)
\] $$
This can be verified by checking
$latex P * A(x) * P^\R{T} = L * L^\R{T}$$.

$head Source$$
$srcfile%example/private/sparse_ad_chol_sp.cpp
	%4%// BEGIN C++%// END C++%1%$$
$end
*/
// BEGIN C++
namespace {
	bool check_sparsity(size_t nx, size_t ny, const CppAD::vector<bool>& s)
	{	bool ok = true;
		//
		// check sparsity for derivative of L_0 (x) = sqrt( x_1 )
		for(size_t j = 0; j < nx; j++)
			ok &= s[ 0 * nx + j ] == (j == 1);
		//
		// check sparsity for derivative of L_1 (x) = sqrt( x_0 )
		for(size_t j = 0; j < nx; j++)
			ok &= s[ 1 * nx + j ] == (j == 0);
		//
		// check sparsity for derivative of L_2 (x) = x_2 / sqrt(x_0)
		for(size_t j = 0; j < nx; j++)
			ok &= s[ 2 * nx + j ] == (j == 0 || j == 2);
		//
		// check for derivative of L_3 (x) = sqrt[ x_3 - x_2^2 / sqrt(x_0) ]
		for(size_t j = 0; j < nx; j++)
			ok &= s[ 3 * nx + j ] == (j == 0 || j == 2 || j == 3 );
		//
		return ok;
	}

}
bool sparse_ad_chol_sp_xam(void)
{	using CppAD::AD;
	typedef CppAD::vector<double>                             d_vector;
	typedef CppAD::vector< AD<double> >                       ad_vector;
	typedef Eigen::SparseMatrix< AD<double>, Eigen::ColMajor> sparse_ad_matrix;
	//
	bool ok     = true;
	// --------------------------------------------------------------------
	size_t nx = 4;
	size_t ny = 4;
	d_vector x(nx), y(ny);
	ad_vector ax(nx), ay(ny);
	ax[0] = x[0] = 2.0;
	ax[1] = x[1] = 1.0;
	ax[2] = x[2] = 0.5;
	ax[3] = x[3] = 3.0;
	// --------------------------------------------------------------------
	// create sparse_ad_cholesky object
	size_t nc = 3;
	sparse_ad_matrix ad_Blow(nc, nc);
	ad_Blow.insert(0,0) = 2.0;
	ad_Blow.insert(1,1) = 0.2;
	ad_Blow.insert(2,0) = 1.0;
	ad_Blow.insert(2,2) = 3.0;
	CppAD::mixed::sparse_ad_cholesky cholesky;
	cholesky.initialize( ad_Blow );
	// ----------------------------------------------------------------------
	// create function object corresponding to L(x)
	CppAD::Independent( ax );
	sparse_ad_matrix ad_Alow(nc, nc);
	ad_Alow.insert(0,0) = ax[0];
	ad_Alow.insert(1,1) = ax[1];
	ad_Alow.insert(2,0) = ax[2];
	ad_Alow.insert(2,2) = ax[3];
	//
	sparse_ad_matrix ad_L;
	cholesky.eval(ad_Alow, ad_L);
	//
	size_t iy = 0;
	for(size_t j = 0; j < nc; j++)
	{	for(sparse_ad_matrix::InnerIterator itr(ad_L, j); itr; ++itr)
			ay[iy++] = itr.value();
	}
	CppAD::ADFun<double> L_fun(ax, ay);
	// ----------------------------------------------------------------------
	ok &= L_fun.Domain() == nx;
	ok &= L_fun.Range()  == ny;
	// ----------------------------------------------------------------------
	// set R to the identity matrix
	size_t q = ny;
	CppAD::vector<bool> r(q * ny);
	for(size_t i = 0; i < q; i++)
	{	for(size_t j = 0; j < ny; j++)
			r[i * ny + j ] = (i == j);
	}
	// compute the sparsity mattern for S(x) = R * L'(x) = L'(x)
	CppAD::vector<bool> s = L_fun.RevSparseJac(q, r);
	ok &= check_sparsity(nx, ny, s);
	// -----------------------------------------------------------------------
	return ok;
}
// END C++
