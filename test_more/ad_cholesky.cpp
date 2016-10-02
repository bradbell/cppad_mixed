// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <Eigen/SparseCore>
# include <Eigen/SparseCholesky>
# include <cppad/cppad.hpp>
# include <cppad/mixed/sparse_mat_info.hpp>
# include <cppad/mixed/sparse_low_tri_sol.hpp>
# include <cppad/mixed/sparse_up_tri_sol.hpp>
# include <cppad/mixed/sparse_scale_diag.hpp>
# include <cppad/mixed/sparse_low2sym.hpp>
# include <cppad/mixed/sparse_mat2low.hpp>
# include <cppad/mixed/sparse_eigen2info.hpp>
# include <iostream>

namespace { // BEGIN_EMPTY_NAMESPACE
// -----------------------------------------------------------------
// define AD, Dynamic, sparse_mat2low, some vector and matrix types
// -----------------------------------------------------------------
using CppAD::AD;
using Eigen::Dynamic;
using Eigen::ColMajor;
using CppAD::mixed::sparse_mat2low;
//
typedef CppAD::vector<size_t>                             s_vector;
typedef CppAD::vector<double>                             d_vector;
typedef CppAD::vector< AD<double> >                       ad_vector;
//
typedef Eigen::PermutationMatrix<Dynamic, Dynamic>        perm_matrix;
typedef Eigen::Matrix<double, Dynamic, Dynamic>           dense_d_matrix;
typedef Eigen::SparseMatrix<double, ColMajor>             sparse_d_matrix;
typedef Eigen::SparseMatrix< AD<double>, ColMajor >       sparse_ad_matrix;
//
/*
sparse_d_matrix sparse_ad2d(const sparse_ad_matrix& amat)
{	sparse_d_matrix result( amat.rows(), amat.cols() );
	for(int k = 0; k < amat.outerSize(); ++k)
	{	for(sparse_ad_matrix::InnerIterator itr(amat, k); itr; ++itr)
			result.insert( itr.row(), itr.col() ) = Value( itr.value() );
	}
	return result;
}
void print_mat(const std::string& label, const sparse_d_matrix& mat)
{	Eigen::Matrix<double, Dynamic, Dynamic> m = mat;
	std::cout << label << "=\n" << m << "\n";
}
*/
// -----------------------------------------------------------------
// get the number of non-zeros corresponding to a matrix
template <class Eigen_sparse_matrix_type>
s_vector get_nnz(const Eigen_sparse_matrix_type& mat)
{	typedef typename Eigen_sparse_matrix_type::InnerIterator iterator;
	size_t nc = size_t( mat.cols() );
	s_vector nnz(nc);
	for(size_t j = 0; j < nc; j++)
	{	nnz[j] = 0;
		for(iterator itr(mat, j); itr; ++itr)
		{	nnz[j]++;
		}
	}
	return nnz;
}
// ------------------------------------------------------------------
// P * M * P^T
sparse_d_matrix pmpt(const perm_matrix& P, const sparse_d_matrix& M)
{	sparse_d_matrix tmp = P * M;
	return tmp * P.transpose();
}
// P^T * M * P
sparse_d_matrix ptmp(const perm_matrix& P, const sparse_d_matrix& M)
{	sparse_d_matrix tmp = P.transpose() * M;
	return tmp * P;
}
// ======================================================================
class atomic_ad_cholesky : public CppAD::atomic_base<double> {
private:
	// -----------------------------------------------------------------
	// data (in order that it is initialized)
	// -----------------------------------------------------------------
	// OK flag
	bool ok_;
	//
	// Number of non-zeros and sparsity pattern for Alow (set by constructor)
	s_vector                      Alow_nnz_;
	CppAD::mixed::sparse_mat_info Alow_pattern_;

	// Object used for Cholesky factorization
	// (analyzePattern by constructor only)
	Eigen::SimplicialLDLT<sparse_d_matrix> ldlt_obj_;

	// Value of the permutation matrix (set by constructor)
	perm_matrix P_;

	// Number of non-zeros and sparsity pattern for L (set by constructor)
	s_vector                      L_nnz_;
	CppAD::mixed::sparse_mat_info L_pattern_;

	// Forward and values for Alow and L (set by forward)
	CppAD::vector<sparse_d_matrix> f_Alow_, f_L_;

public:
	// -----------------------------------------------------------------
	// normal member functions
	// -----------------------------------------------------------------
	// status flag
	bool ok(void)
	{	return ok_; }
	// -----------------------------------------------------------------
	// constructor
	atomic_ad_cholesky(const sparse_d_matrix& Alow )
	: CppAD::atomic_base<double> (
		"atomic_ad_cholesky",
		CppAD::atomic_base<double>::set_sparsity_enum
	)
	{	ok_ = true;
		//
		// Step 1: Set Alow_nnz_ and Alow_pattern_
		Alow_nnz_     = get_nnz(Alow);
		Alow_pattern_ = CppAD::mixed::sparse_eigen2info(Alow);
		//
		// Step 2: analyze the sparsity pattern
		ldlt_obj_.analyzePattern( Alow );
		//
		// Step 3: Compute the Cholesky factor for this Alow
		// and the permutation all Alow (used with this object).
		ldlt_obj_.factorize( Alow );
		ok_ &= ldlt_obj_.info() == Eigen::Success;
		//
		// Step 4: Retrieve the permutation P_;
		P_   = ldlt_obj_.permutationP();
		//
		// Step 5: Set L_nnz_ and L_pattern_
		sparse_d_matrix L  =  ldlt_obj_.matrixL();
		L_nnz_             = get_nnz(L);
		L_pattern_         = CppAD::mixed::sparse_eigen2info(L);
	}
	// -----------------------------------------------------------------
	// user AD version of atomic Cholesky factorization
	void ad(
		const sparse_ad_matrix& aAlow  ,
		sparse_ad_matrix&       aL     )
	{	// -----------------------------------------------------------
		// packed version of Alow
		size_t nc = Alow_nnz_.size();
		size_t nx = Alow_pattern_.row.size();
		assert( nc == size_t( aAlow.rows() ) );
		assert( nc == size_t( aAlow.cols() ) );
		ad_vector ax( nx );
		size_t index = 0;
		for(size_t j = 0; j < nc; j++)
		{	for(sparse_ad_matrix::InnerIterator itr(aAlow, j); itr; ++itr)
			{	assert( Alow_pattern_.row[index] == size_t( itr.row() ) );
				assert( Alow_pattern_.col[index] == size_t( itr.col() ) );
				ax[ index ] = itr.value();
				++index;
			}
		}
		assert( index == nx );
		// -------------------------------------------------------------------
		// make call to packed vector verison of the atomic function
		size_t ny = L_pattern_.row.size();
		ad_vector ay( ny );
		(*this)(ax, ay);
		// -------------------------------------------------------------------
		// unpack ay into aL
		aL.resize(nc, nc);
		index = 0;
		for(size_t ell = 0; ell < ny; ell++)
		{	size_t i = L_pattern_.row[ell];
			size_t j = L_pattern_.col[ell];
			aL.insert(i, j) = ay[ index++ ];
		}
		assert( index == ny );
		return;
	}
private:
	// -----------------------------------------------------------------
	// virtual member functions
	// ------------------------------------------------------------------
	// CppAD forward mode for this operation
	virtual bool forward(
		// lowest order Taylor coefficient we are evaluating
		size_t                          p ,
		// highest order Taylor coefficient we are evaluating
		size_t                          q ,
		// which components of x are variables
		const CppAD::vector<bool>&      vx ,
		// which components of y are variables
		CppAD::vector<bool>&            vy ,
		// tx [ j * (q+1) + k ] is x_j^k
		const CppAD::vector<double>&    tx ,
		// ty [ i * (q+1) + k ] is y_i^k
		CppAD::vector<double>&          ty )
	{	typedef typename sparse_d_matrix::InnerIterator iterator;
		//
		assert( p <= q );
		//
		size_t n_order = q + 1;
		size_t nx      = Alow_pattern_.row.size();
		size_t ny      = L_pattern_.row.size();
		size_t nc      = Alow_nnz_.size();
		assert( vx.size() == 0 || nx == vx.size() );
		assert( vx.size() == 0 || ny == vy.size() );
		assert( nx * n_order == tx.size() );
		assert( ny * n_order == ty.size() );
		// -------------------------------------------------------------------
		// make sure f_Alow_ and f_L_ are large enough
		assert( f_Alow_.size() == f_L_.size() );
		if( f_Alow_.size() < n_order )
		{	f_Alow_.resize(n_order);
			f_L_.resize(n_order);
		}
		// -------------------------------------------------------------------
		// unpack tx into f_Alow_
		for(size_t k = 0; k < n_order; k++)
		{	// unpack Alow values for this order
			f_Alow_[k].resize(nc, nc);
			f_Alow_[k].reserve(Alow_nnz_);
			for(size_t ell = 0; ell < Alow_pattern_.row.size(); ell++)
			{	size_t i = Alow_pattern_.row[ell];
				size_t j = Alow_pattern_.col[ell];
				f_Alow_[k].insert(i, j) = tx[ ell * n_order + k ];
			}
		}
		// -------------------------------------------------------------------
		// for orders less than p, unpack ty into f_L_
		for(size_t k = 0; k < p; k++)
		{	// unpack f_L_ values for this order
			assert( L_pattern_.row.size() == ny );
			f_L_[k].resize(nc, nc);
			f_L_[k].reserve(L_nnz_);
			for(size_t index = 0; index < ny; index++)
			{	size_t i = L_pattern_.row[index];
				size_t j = L_pattern_.col[index];
				f_L_[k].insert(i, j) = ty[ index * n_order + k];
			}
		}
		if( p == 0 )
		{	// compute result for zero order
			ldlt_obj_.factorize( f_Alow_[0] );
			if( ldlt_obj_.info() != Eigen::Success )
			{	ok_ = false;
				return false;
			}
			Eigen::VectorXd D2 = sqrt( ldlt_obj_.vectorD().array() ).matrix();
			f_L_[0]            =  ldlt_obj_.matrixL() * D2.asDiagonal();
			assert( P_.indices() == ldlt_obj_.permutationP().indices() );
		}
		Eigen::SparseMatrix<double, Eigen::RowMajor> L0 = f_L_[0];
		// compute result for orders greater than or equal 1 and p.
		for(size_t k = std::max(p, size_t(1)); k < n_order; k++)
		{	// convert Alow to A_k
			sparse_d_matrix A_k = CppAD::mixed::sparse_low2sym(f_Alow_[k]);
			//
			// initialize sum as E_k =  P * A_k * P.transpose()
			sparse_d_matrix f_sum  = pmpt(P_, A_k);
			//
			// compute E_k - B_k
			for(size_t ell = 1; ell < k; ell++)
				f_sum -= f_L_[ell] * f_L_[k-ell].transpose();
			// compute L_0^{-1} * (E_k - B_k)
			sparse_d_matrix tmp1 = CppAD::mixed::sparse_low_tri_sol(L0, f_sum);
			//
			// compute L_0^{-1} * (E_k - B_k) * L_0^{-T}
			sparse_d_matrix tmp2 = CppAD::mixed::sparse_low_tri_sol(
				L0, tmp1.transpose()
			).transpose();
			//
			// divide the diagonal by 2
			CppAD::mixed::sparse_scale_diag(0.5, tmp2);
			//
			// low[ L_0^{-1} * (E_k - B_k) * L_0^{-T} ]
			// L_k = L_0 * low[  L_0^{-1} * (E_k - B_k) * L_0^{-T} ]
			f_L_[k] = L0 * sparse_mat2low(tmp2);
		}
		// -------------------------------------------------------------------
		// pack f_L_ into ty
		// -------------------------------------------------------------------
		for(size_t k = p; k < n_order; k++)
		{	size_t index = 0;
			// repacking is harder because some of the computed terms
			// may be zero and may not appear in f_L_[k];
			for(size_t j = 0; j < nc; j++)
			{	//
				// initialize this column as zero
				for(size_t ell = 0; ell < L_nnz_[j]; ell++)
					ty[ (index + ell) * n_order + k ] = 0.0;
				//
				// fill in non-zero entries
				size_t ell = 0;
				for(iterator itr(f_L_[k], j); itr; ++itr)
				{	size_t i    = L_pattern_.row[index + ell];
					while( i < size_t(itr.row()) )
					{	++ell;
						i = L_pattern_.row[index + ell];
					}
					assert( i == size_t( itr.row() ) );
					assert( L_pattern_.col[index + ell] == j );
					ty[ (index + ell) * n_order + k ] = itr.value();
				}
				index += L_nnz_[j];
			}
		}
		// -------------------------------------------------------------------
		// check if we are computing vy
		if( vx.size() == 0 )
			return true;
		// -------------------------------------------------------------------
		// This is a very dumb algorithm that over estimates which elements
		// of L are variables. 2DO: create a much better estimate
		bool var = false;
		for(size_t j = 0; j < nx; j++)
			var |= vx[j];
		for(size_t i = 0; i < ny; i++)
			vy[i] = var;
		//
		return true;
	}
	// ------------------------------------------------------------------
	// reverse mode routine called by CppAD
	virtual bool reverse(
		// highest order Taylor coefficient that we are computing derivative of
		size_t                     q ,
		// forward mode Taylor coefficients for x variables
		const CppAD::vector<double>&     tx ,
		// forward mode Taylor coefficients for y variables
		const CppAD::vector<double>&     ty ,
		// upon return, derivative of G[ F[ {x_j^k} ] ] w.r.t {x_j^k}
		CppAD::vector<double>&           px ,
		// derivative of G[ {y_i^k} ] w.r.t. {y_i^k}
		const CppAD::vector<double>&     py
	)
	{	typedef typename sparse_d_matrix::InnerIterator iterator;
		size_t n_order = q + 1;
		size_t nx      = Alow_pattern_.row.size();
		size_t ny      = L_pattern_.row.size();
		size_t nc      = Alow_nnz_.size();
		//
		assert( nx * n_order == tx.size() );
		assert( ny * n_order == ty.size() );
		assert( px.size()    == tx.size() );
		assert( py.size()    == ty.size() );
		//
		assert( f_Alow_.size() == f_L_.size() );
		assert( f_Alow_.size() >= n_order );
		// -------------------------------------------------------------------
		// declare r_Alow, r_L
		CppAD::vector<sparse_d_matrix> r_Alow(n_order), r_L(n_order);
		// -------------------------------------------------------------------
		// unpack tx into f_Alow_
		assert( Alow_pattern_.row.size() == nx );
		for(size_t k = 0; k < n_order; k++)
		{	// unpack Alow_ for this order
			f_Alow_[k].resize(nc, nc);
			f_Alow_[k].reserve(Alow_nnz_);
			for(size_t index = 0; index < nx; index++)
			{	size_t i = Alow_pattern_.row[index];
				size_t j = Alow_pattern_.col[index];
				f_Alow_[k].insert(i, j) = tx [ index * n_order + k ];
			}
		}
		// -------------------------------------------------------------------
		// for orders less than p, unpack ty into f_L_
		assert( L_pattern_.row.size() == ny );
		for(size_t k = 0; k < n_order; k++)
		{	// unpack f_L_ values for this order
			f_L_[k].setZero();
			f_L_[k].reserve(L_nnz_);
			for(size_t index = 0; index < ny; index++)
			{	size_t i = L_pattern_.row[index];
				size_t j = L_pattern_.col[index];
				f_L_[k].insert(i, j) = ty[ index * n_order + k];
			}
		}
		// -------------------------------------------------------------------
		// unpack py into r_L
		assert( L_pattern_.row.size() == ny );
		for(size_t k = 0; k < n_order; k++)
		{	r_L[k].resize(nc, nc);
			r_L[k].reserve(L_nnz_);
			for(size_t ell = 0; ell < ny; ell++)
			{	size_t i = L_pattern_.row[ell];
				size_t j = L_pattern_.col[ell];
				r_L[k].insert(i, j) = py[ ell * n_order + k];
			}
		}
		// -------------------------------------------------------------------
		// initialize r_Alow as zero
		for(size_t k = 0; k < n_order; k++)
			r_Alow[k].resize(nc, nc);
		// -------------------------------------------------------------------
		// Cholesky factorization
		Eigen::SparseMatrix<double, Eigen::RowMajor> L0 = f_L_[0];
		//
		// start at highest order and go down
		for(size_t k1 = n_order; k1 > 1; k1--)
		{	size_t k = k1 - 1;
			// L_0^T * bar{L}_k
			sparse_d_matrix tmp1 = L0.transpose() * r_L[k];
			//
			// low[ L_0^T * bar{L}_k ]
			CppAD::mixed::sparse_scale_diag(0.5, tmp1);
			sparse_d_matrix tmp2 = sparse_mat2low(tmp1);
			//
			// L_0^{-T} * low[ L_0^T * bar{L}_k ]
			tmp1 = CppAD::mixed::sparse_up_tri_sol(L0.transpose(), tmp2);
			//
			// L_0^{-T} * low[ L_0^T * bar{L}_k ]^T * L_0^{-1}
			sparse_d_matrix Mk = CppAD::mixed::sparse_up_tri_sol(
				L0.transpose(), tmp1.transpose()
			).transpose();
			//
			// remove Lk, \bar{Alow}_k += P^T * M0 * P
			sparse_d_matrix barB_k = CppAD::mixed::sparse_low2sym(Mk);
			r_Alow[k]             += sparse_mat2low( ptmp(P_, barB_k) );
			// compute barB_k
			barB_k                 = -1.0 * barB_k;
			//
			// remove C_k using
			// 2 * lower[ bar{B}_k L_k ]
			tmp1    = barB_k * f_L_[k];
			r_L[0] += 2.0 * sparse_mat2low( tmp1 );
			//
			// remove B_k
			for(size_t ell = 1; ell < k; ell++)
			{	// bar{L}_ell = 2 * lower( bar{B}_k * L_{k-ell} )
				tmp1      = barB_k * f_L_[k-ell];
				r_L[ell] += 2.0 * sparse_mat2low( tmp1 );
			}
		}
		// L_0^T * bar{L}_0
		sparse_d_matrix tmp1 = L0.transpose() * r_L[0];
		//
		// low[ L_0^T * bar{L}_0 ]
		CppAD::mixed::sparse_scale_diag(0.5, tmp1);
		sparse_d_matrix tmp2 = sparse_mat2low( tmp1 );
		//
		// L_0^{-T} low[ L_0^T * bar{L}_0 ]
		tmp1 = CppAD::mixed::sparse_up_tri_sol(L0.transpose(), tmp2);
		//
		// M_0 = L_0^{-T} low[ L_0^T * bar{L}_0 ]^T L_0^{-1}
		sparse_d_matrix M0 = CppAD::mixed::sparse_up_tri_sol(
			L0.transpose(), tmp1.transpose()
		);
		// remove L0, \bar{Alow}_0 += 2.0 * low[ P^T * M0 * P ]
		tmp1       = 2.0 * ptmp(P_, M0);
		CppAD::mixed::sparse_scale_diag(0.5, tmp1);
		r_Alow[0] += sparse_mat2low( tmp1 );
		// ------------------------------------------------------------------
		// pack r_Alow into px
		for(size_t k = 0; k < n_order; k++)
		{	size_t index = 0;
			// repacking is harder because some of the computed terms
			// may be zero and not appear in r_Alow[k]
			for(size_t j = 0; j < nc; j++)
			{	// initialize this column as zero
				for(size_t ell = 0; ell < Alow_nnz_[j]; ell++)
					px[ (index + ell) * n_order + k ] = 0.0;

				// fill in the non-zero entries
				size_t ell = 0;
				for(iterator itr(r_Alow[k], j); itr; ++itr)
				{	size_t i = Alow_pattern_.row[index + ell];
					while( i < size_t( itr.row() ) )
					{	++ell;
						i = Alow_pattern_.row[index + ell];
					}
					assert( i == size_t( itr.row() ) );
					assert( Alow_pattern_.col[index + ell] == j );
					px[ (index + ell) * n_order + k ] = itr.value();
				}
				index += Alow_nnz_[j];
			}
		}
		//
		return true;
	}
}; // END_ATOMIC_AD_CHOLESKY

} // END_EMPTY_NAMESPACE

/*
$begin ad_cholesky.cpp$$

$section Sparse Atomic AD Cholesky Factorization: Example and Test$$

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
The leading princial minors of this matrix are
$latex x_0$$,
$latex x_0 x_1$$,
$latex x_0 ( x_1 x_2 - x_1 x_1 )$$,
If all these minors are positive, the matrix $latex A(x)$$ is
positive definite.

$end
*/

bool ad_cholesky(void)
{	bool ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();
	// --------------------------------------------------------------------
	size_t nx = 3;
	size_t ny = 1;
	d_vector x(nx), y(ny);
	x[0] = 1.0;
	x[1] = 0.5;
	x[2] = 2.0;
	// --------------------------------------------------------------------
	// create ad_cholesky object
	size_t nc = 3;
	sparse_d_matrix Blow(nc, nc);
	Blow.insert(0,0) = x[0];
	Blow.insert(2,0) = x[1];
	Blow.insert(1,1) = x[1];
	Blow.insert(2,2) = x[2];
	atomic_ad_cholesky cholesky( Blow );
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
	sparse_ad_matrix aAlow(nc, nc);
	aAlow.insert(0,0) = ax[0];
	aAlow.insert(2,0) = ax[1];
	aAlow.insert(1,1) = ax[1];
	aAlow.insert(2,2) = ax[2];
	//
	// compute the Choleksy factorization of A
	sparse_ad_matrix aL;
	cholesky.ad(aAlow, aL);
	//
	// diagonal of L
	Eigen::Matrix< AD<double> , Dynamic , 1 > D = aL.diagonal();
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
