// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-20 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/typedef.hpp>
# include "../sparse_ad_cholesky.hpp"

/*
$begin sparse_ad_chol_sp2.cpp$$
$spell
	Cholesky
$$

$section Sparse AD Cholesky Sparsity Calculation: Example and Test$$

$head Problem$$
We are given the function
$latex A : \B{R}^6 \rightarrow \B{R}^{4 \times 4}$$
defined by
$latex \[
	A(x) = \left( \begin{array}{cccc}
		x_0 & x_1  & 0   & 0   \\
		x_1 & x_2  & x_3 & 0   \\
		0   & x_3  & x_4 & x_5 \\
		0   & 0    & x_5 & x_6
	\end{array} \right)
\] $$
Note that the subscripts are in column major order for the
lower triangle of $latex A(x)$$.

$head Permutation$$
The fill reducing permutation
$cref/P/sparse_ad_cholesky/Notation/P/$$
transposes indices zero and one; i.e.
$latex \[
	P= \left( \begin{array}{cccc}
		1   & 0    & 0  & 0  \\
		0   & 0    & 0  & 1  \\
		0   & 1    & 0  & 0  \\
		0   & 0    & 1  & 0  \\
	\end{array} \right)
	\; , \;
	P A(x) = \left( \begin{array}{cccc}
		x_0 & x_1  & 0   & 0   \\
		0   & 0    & x_5 & x_6 \\
		x_1 & x_2  & x_3 & 0   \\
		0   & x_3  & x_4 & x_5
	\end{array} \right)
	\; , \;
	P A(x) P^T = \left( \begin{array}{cccc}
		x_0 & 0   & x_1 & 0   \\
		0   & x_6 & 0   & x_5 \\
		x_1 & 0   & x_2 & x_3 \\
		0   & x_5 & x_3 & x_4
	\end{array} \right)
\] $$

$head Cholesky Factor$$
Define the function
$latex \[
	c(x) = x_3 / \sqrt{x_2 - x_1^2/ x_0}
\] $$
The Cholesky factor
$cref/L/sparse_ad_cholesky/Notation/L/$$ is
$latex \[
L(x) = \left(
\begin{array}{cccc}
\sqrt{x_0}       & 0                &  0                      & 0 \\
0                & \sqrt{x_6}       &  0                      & 0 \\
x_1 / \sqrt{x_0} & 0                & \sqrt{x_2 - x_1^2/ x_0} & 0 \\
0                & x_5 / \sqrt{x_6} & c(x) & \sqrt{x_4 - x_5^2 / x_6 - c(x)^2 }
\end{array}
\right)
\] $$
This can be verified by checking
$latex P * A(x) * P^\R{T} = L * L^\R{T}$$.

$head Source$$
$srcthisfile%4%// BEGIN C++%// END C++%1%$$
$end
*/

// include tests that are not yet working
# define CPPAD_MIXED_INCLUDE_TEST_NOT_YET_WORKING 0

// BEGIN C++
namespace {
	bool check_fun(
		const CppAD::vector<double>& x ,
		const CppAD::vector<double>& y )
	{	bool ok = true;
		double eps = 100. * std::numeric_limits<double>::epsilon();
		using CppAD::NearEqual;
		using std::sqrt;
		//
		ok &= y.size() == 7;
		//
		double check = sqrt( x[0] );
		ok &= NearEqual(y[0], check, eps, eps);
		//
		check = x[1] / sqrt( x[0] );
		ok &= NearEqual(y[1], check, eps, eps);
		//
		check = sqrt( x[6] );
		ok &= NearEqual(y[2], check, eps, eps);
		//
		check = x[5] / sqrt( x[6] );
		ok &= NearEqual(y[3], check, eps, eps);
		//
		check = sqrt( x[2] - x[1]* x[1] / x[0] );
		ok &= NearEqual(y[4], check, eps, eps);
		//
		double c = x[3] / sqrt( x[2] - x[1]* x[1] / x[0] );
		ok &= NearEqual(y[5], c, eps, eps);
		//
		check = sqrt( x[4] - x[5] * x[5] / x[6] - c * c );
		ok &= NearEqual(y[6], check, eps, eps);
		//
		return ok;
	}
	bool check_jac_sparsity(size_t nx, const CppAD::vector<bool>& s)
	{	bool ok = true;
		//
		// derivative of L_0 (x) = sqrt( x_0 )
		for(size_t j = 0; j < nx; j++)
			ok &= s[ 0 * nx + j ] == (j == 0);
		//
		// derivative of L_1 (x) = x_1 / sqrt(x_0)
		for(size_t j = 0; j < nx; j++)
			ok &= s[ 1 * nx + j ] == (j == 0 || j == 1) ;
		//
		// derivative of L_2 (x) = sqrt(x_6)
		for(size_t j = 0; j < nx; j++)
			ok &= s[ 2 * nx + j ] == (j == 6);
		//
		// derivative of L_3 (x) = x_5 / sqrt(x_6)
		for(size_t j = 0; j < nx; j++)
			ok &= s[ 3 * nx + j ] == (j == 5 || j == 6 );
		//
		// derivative of L_4 (x) = sqrt[ x_2 - x_1^2 / x_0 ]
		for(size_t j = 0; j < nx; j++)
			ok &= s[ 4 * nx + j ] == ( j == 0 || j == 1 || j == 2 );
		//
		// derivative of L_5 (x) = c(x) = x_3 / sqrt[ x_2 - x_1^2 / x_0 ]
		for(size_t j = 0; j < nx; j++)
			ok &= s[ 5 * nx + j ] == ( j <= 3 );
		//
		// derivative of L_6 (x) = sqrt[ x_4 - x_5^2 / x_6  - c(x)^2 ]
		for(size_t j = 0; j < nx; j++)
			ok &= s[ 6 * nx + j ] == (j <= 6);
		//
		return ok;
	}
	bool check_hes_sparsity(size_t nx, size_t k, const CppAD::vector<bool>& h)
	{	bool ok = true;
		switch(k)
		{	case 0:
			// L_0 (x) = sqrt( x_0 )
			for(size_t i = 0; i < nx; i++)
			{	for(size_t j = 0; j < nx; j++)
				{	bool check = (i == 0) && (j == 0);
					ok &= h[ i * nx + j] == check;
				}
			}
			break;

			case 1:
			// L_1 (x) = x_1 / sqrt(x_0)
			for(size_t i = 0; i < nx; i++)
			{	for(size_t j = 0; j < nx; j++)
				{	bool check = false;
					check |= (i == 0 && j == 0);
					check |= (i == 0 && j == 1) || (i == 1 && j == 0);
					ok &= h[ i * nx + j] == check;
				}
			}
			break;

			case 2:
			// L_2 (x) = sqrt(x_6)
			for(size_t i = 0; i < nx; i++)
			{	for(size_t j = 0; j < nx; j++)
				{	bool check = false;
					check |= (i == 6 && j == 6);
					ok &= h[ i * nx + j] == check;
				}
			}
			break;

			case 3:
			// L_3 (x) = x_5 / sqrt(x_6)
			for(size_t i = 0; i < nx; i++)
			{	for(size_t j = 0; j < nx; j++)
				{	bool check = false;
					check |= (i == 6 && j == 6);
					check |= (i == 5 && j == 6);
					check |= (j == 5 && i == 6);
					ok &= h[ i * nx + j] == check;
				}
			}
			break;

			case 4:
			// L_4 (x) = sqrt[ x_2 - x_1^2 / x_0 ]
			for(size_t i = 0; i < nx; i++)
			{	for(size_t j = 0; j < nx; j++)
				{	bool check = false;
					check |= (i <= 2 && j <= 2);
					ok &= h[ i * nx + j] == check;
				}
			}
			break;

			case 5:
			// L_5 (x) = x_3 / sqrt[ x_2 - x_1^2 / x_0 ]
			for(size_t i = 0; i < nx; i++)
			{	for(size_t j = 0; j < nx; j++)
				{	bool check = false;
					check |= (i <= 2 && j <= 2);
					check |= (i <= 2 && j == 3);
					check |= (i == 3 && j <= 2);
					ok &= h[ i * nx + j] == check;
				}
			}
			break;

			case 6:
			// L_6 (x) = sqrt[ x_4 - x_5^2 / x_6  - c(x)^2 ]
			for(size_t i = 0; i < nx; i++)
			{	for(size_t j = 0; j < nx; j++)
				{	bool check = false;
					check |= (i <= 6 && j <= 6);
					ok &= h[ i * nx + j] == check;
				}
			}
			break;

			default:
			ok = false;
		}
		return ok;
	}
}
bool sparse_ad_chol_sp2(void)
{	using CppAD::AD;
	using CppAD::mixed::d_vector;
	typedef CppAD::vector< AD<double> >                       ad_vector;
	typedef Eigen::SparseMatrix< AD<double>, Eigen::ColMajor> sparse_ad_matrix;
	using Eigen::Dynamic;
	//
	bool ok     = true;
	// --------------------------------------------------------------------
	size_t nx = 7;
	size_t ny = 7;
	d_vector x(nx), y(ny);
	ad_vector ax(nx), ay(ny);
	ax[0] = x[0] = 1.0;
	ax[1] = x[1] = 0.1;
	ax[2] = x[2] = 2.0;
	ax[3] = x[3] = 0.3;
	ax[4] = x[4] = 4.0;
	ax[5] = x[5] = 0.5;
	ax[6] = x[6] = 6.0;
	// --------------------------------------------------------------------
	// create sparse_ad_cholesky object
	size_t nc = 4;
	sparse_ad_matrix ad_Blow;
	ad_Blow.resize( int(nc), int(nc) );
	ad_Blow.insert(0,0) = x[0];
	ad_Blow.insert(1,0) = x[1];
	ad_Blow.insert(1,1) = x[2];
	ad_Blow.insert(2,1) = x[3];
	ad_Blow.insert(2,2) = x[4];
	ad_Blow.insert(3,2) = x[5];
	ad_Blow.insert(3,3) = x[6];
	CppAD::mixed::sparse_ad_cholesky cholesky;
	cholesky.initialize( ad_Blow );
	// ----------------------------------------------------------------------
	// check that the permutation matrix P is as expected
	Eigen::PermutationMatrix<Dynamic, Dynamic> P = cholesky.permutation();
	Eigen::Matrix<size_t, Dynamic, 1>  indices(nc);
	for(size_t i = 0; i < nc; i++)
			indices[i] = i;
	indices = P * indices;
	ok &= indices[0] == 0;
	ok &= indices[1] == 3;
	ok &= indices[2] == 1;
	ok &= indices[3] == 2;
	// ----------------------------------------------------------------------
	// create function object corresponding to L(x)
	CppAD::Independent( ax );
	sparse_ad_matrix ad_Alow;
	ad_Alow.resize( int(nc), int(nc) );
	ad_Alow.insert(0,0) = ax[0];
	ad_Alow.insert(1,0) = ax[1];
	ad_Alow.insert(1,1) = ax[2];
	ad_Alow.insert(2,1) = ax[3];
	ad_Alow.insert(2,2) = ax[4];
	ad_Alow.insert(3,2) = ax[5];
	ad_Alow.insert(3,3) = ax[6];
	//
	sparse_ad_matrix ad_L;
	cholesky.eval(ad_Alow, ad_L);
	//
	size_t iy = 0;
	for(size_t j = 0; j < nc; j++)
	{	for(sparse_ad_matrix::InnerIterator itr(ad_L, int(j)); itr; ++itr)
			ay[iy++] = itr.value();
	}
	CppAD::ADFun<double> L_fun(ax, ay);
	ok &= L_fun.Domain() == nx;
	ok &= L_fun.Range()  == ny;
	// ----------------------------------------------------------------------
	y   = L_fun.Forward(0, x);
	ok &= check_fun(x, y);
	// ----------------------------------------------------------------------
	// Check rev_sparse_jac
	//
	// set R to the identity matrix
	size_t q = ny;
	CppAD::vector<bool> r(q * ny);
	for(size_t i = 0; i < q; i++)
	{	for(size_t j = 0; j < ny; j++)
			r[i * ny + j ] = (i == j);
	}
	// compute the sparsity pattern for S(x) = R * L'(x) = L'(x)
	CppAD::vector<bool> s = L_fun.RevSparseJac(q, r);
	ok &= check_jac_sparsity(nx, s);
	// ----------------------------------------------------------------------
	// Check for_sparse_jac
	//
	// set R to the identity matrix
	q = nx;
	r.resize(nx * q);
	for(size_t i = 0; i < nx; i++)
	{	for(size_t j = 0; j < q; j++)
			r[i * q + j ] = (i == j);
	}
	// compute the sparsity pattern for S(x) = L'(x) * R = L'(x)
	s = L_fun.ForSparseJac(q, r);
	ok &= check_jac_sparsity(nx, s);
	// -----------------------------------------------------------------------
	// Check rev_sparse_hes
	s.resize(ny);
	for(size_t i = 0; i < ny; i++)
	{	// set s
		for(size_t j = 0; j < ny; j++)
			s[j] = (i == j);
		//
		// compute sparsity pattern for Hessian of s * L(x)
		CppAD::vector<bool> h = L_fun.RevSparseHes(q, s);
		//
		// check sparsity for this component function
		ok &= check_hes_sparsity(nx, i, h);
	}
	return ok;
}
// END C++
