/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/typedef.hpp>
# include "../sparse_ad_cholesky.hpp"

/*
$begin sparse_ad_chol_eval.cpp$$
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
$latex x_0 x_1 x_2 - x_1 x_1 x_1$$,
If all these minors are positive, the matrix $latex A(x)$$ is
positive definite.

$head Source$$
$srcfile%cholesky/example/sparse_ad_chol_eval.cpp
	%4%// BEGIN C++%// END C++%1%$$
$end
*/
// BEGIN C++
bool sparse_ad_chol_eval(void)
{	using CppAD::AD;
	using CppAD::mixed::d_vector;
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
	Eigen::SparseMatrix< AD<double>, Eigen::ColMajor> ad_Blow;
	ad_Blow.resize(int(nc), int(nc));
	ad_Blow.insert(0,0) = x[0];
	ad_Blow.insert(2,0) = x[1];
	ad_Blow.insert(1,1) = x[1];
	ad_Blow.insert(2,2) = x[2];
	CppAD::mixed::sparse_ad_cholesky cholesky;
	cholesky.initialize( ad_Blow );
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
	Eigen::SparseMatrix< AD<double>, Eigen::ColMajor> ad_Alow;
	ad_Alow.resize(int(nc), int(nc));
	ad_Alow.insert(0,0) = ax[0];
	ad_Alow.insert(2,0) = ax[1];
	ad_Alow.insert(1,1) = ax[1];
	ad_Alow.insert(2,2) = ax[2];
	//
	// compute the Choleksy factorization of A
	Eigen::SparseMatrix< AD<double>, Eigen::ColMajor> ad_L;
	cholesky.eval(ad_Alow, ad_L);
	ok &= size_t(ad_L.rows()) == nc;
	ok &= size_t(ad_L.cols()) == nc;
	//
	// diagonal of L
	Eigen::Matrix< AD<double> , Eigen::Dynamic , 1 > ad_D = ad_L.diagonal();
	//
	// product of diagonal elements of L
	AD<double> ad_p = 1.0;
	for(size_t j = 0; j < nc; j++)
		ad_p *= ad_D[j];
	//
	// calculate the determinant of A using the Cholesky factorization
	ad_vector ay(ny);
	ay[0] = ad_p * ad_p;
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
