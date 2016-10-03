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
$begin sparse_ad_cholesky_ad_xam.cpp$$
$spell
	Cholesky
$$

$section Sparse AD Cholesky Factorization: Example and Test$$

$head Problem$$
We are given the function $latex A : \B{R}^3 \rightarrow \B{R}^{3 \times 3}$$
defined by
$latex \[
	A(x) = \left( \begin{array}{ccc}
		x_0 & 0    & x_1  \\
		0   & x_1  & 0   \\
		x_1 & 0    & x_2
	\end{array} \right)
\] $$
The leading principal minors of this matrix are
$latex x_0$$,
$latex x_0 x_1$$,
$latex x_0 ( x_1 x_2 - x_1 x_1 )$$,
If all these minors are positive, the matrix $latex A(x)$$ is
positive definite.

$head Source$$
$srcfile%example/private/sparse_ad_cholesky_ad.cpp
	%4%// BEGIN C++%// END C++%1%$$
$end
*/
// BEGIN C++
bool sparse_ad_cholesky_ad_xam(void)
{	using CppAD::AD;
	typedef CppAD::vector<double>         d_vector;
	typedef CppAD::vector< AD<double> >   ad_vector;
	//
	bool ok     = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();
	// --------------------------------------------------------------------
	size_t nx = 3;
	size_t ny = 1;
	d_vector x(nx), y(ny);
	x[0] = 1.0;
	x[1] = 0.5;
	x[2] = 2.0;
	// --------------------------------------------------------------------
	// create sparse_ad_cholesky object
	size_t nc = 3;
	Eigen::SparseMatrix<double, Eigen::ColMajor> Blow(nc, nc);
	Blow.insert(0,0) = x[0];
	Blow.insert(2,0) = x[1];
	Blow.insert(1,1) = x[1];
	Blow.insert(2,2) = x[2];
	CppAD::mixed::sparse_ad_cholesky cholesky( Blow );
	//
	// ----------------------------------------------------------------------
	// Create function object corresponding to f(x)
	//
	// Independent variables
	ad_vector ax(nx);
	ax[0] = x[0];
	ax[1] = x[1];
	ax[2] = x[2];
	CppAD::Independent( ax );
	//
	// Lower triangle of symmetric matrix with same sparsity pattern as B
	Eigen::SparseMatrix< AD<double>, Eigen::ColMajor> aAlow(nc, nc);
	aAlow.insert(0,0) = ax[0];
	aAlow.insert(2,0) = ax[1];
	aAlow.insert(1,1) = ax[1];
	aAlow.insert(2,2) = ax[2];
	//
	// compute the Choleksy factorization of A
	Eigen::SparseMatrix< AD<double>, Eigen::ColMajor> aL(nc, nc);
	cholesky.ad(aAlow, aL);
	//
	// diagonal of L
	Eigen::Matrix< AD<double> , Eigen::Dynamic , 1 > D = aL.diagonal();
	//
	// product of diagonal elements of L
	AD<double> p = 1.0;
	for(size_t j = 0; j < nc; j++)
		p *= D[j];
	//
	// calculate the determinant of A using the Cholesky factorization
	ad_vector ay(ny);
	ay[0] = p * p;
	//
	// f(x) = det[ A(x) ]
	CppAD::ADFun<double> f(ax, ay);
	// ----------------------------------------------------------------------
	// Test zero order forward
	y    = f.Forward(0, x);
	double check = x[0] * x[1] * x[2] - x[1] * x[1] * x[1];
	ok          &= CppAD::NearEqual(y[0], check, eps, eps );
	// -----------------------------------------------------------------------
	// Test first order forward
	d_vector x1(nx), y1(ny);
	//
	// partial w.r.t. x[0]
	x1[0]  = 1.0;
	x1[1]  = 0.0;
	x1[2]  = 0.0;
	y1     = f.Forward(1, x1);
	double f_x0 = x[1] * x[2];
	ok        &= CppAD::NearEqual(y1[0], f_x0, eps, eps);
	//
	// partial w.r.t. x[2]
	x1[0]  = 0.0;
	x1[2]  = 1.0;
	y1     = f.Forward(1, x1);
	double f_x2 = x[0] * x[1];
	ok        &= CppAD::NearEqual(y1[0], f_x2, eps, eps);
	//
	// partial w.r.t. x[1]
	x1[2]  = 0.0;
	x1[1]  = 1.0;
	y1     = f.Forward(1, x1);
	double f_x1 = x[0] *  x[2] - 3.0 * x[1] * x[1];
	ok    &= CppAD::NearEqual(y1[0], f_x1, eps, eps);
	// -----------------------------------------------------------------------
	// Test second order forward
	d_vector x2(nx), y2(ny);
	//
	// second partial w.r.t x[1]
	x2[0]  = 0.0;
	x2[1]  = 0.0;
	x2[2]  = 0.0;
	y2     = f.Forward(2, x2);
	double f_x11  = - 6.0 * x[1];
	ok   &= CppAD::NearEqual(y2[0], f_x11 / 2.0, eps, eps);
	// -----------------------------------------------------------------------
	// Test first order reverse
	d_vector w(ny), d1w(nx);
	w[0] = 1.0;
	d1w  = f.Reverse(1, w);
	ok  &= CppAD::NearEqual(d1w[0], f_x0, eps, eps);
	ok  &= CppAD::NearEqual(d1w[1], f_x1, eps, eps);
	ok  &= CppAD::NearEqual(d1w[2], f_x2, eps, eps);
	// -----------------------------------------------------------------------
	// check for any errors during use of cholesky
	ok &= cholesky.ok();
	// -----------------------------------------------------------------------
	return ok;
}
// END C++
