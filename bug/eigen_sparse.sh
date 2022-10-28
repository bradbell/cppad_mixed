# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-22 Bradley M. Bell
# ----------------------------------------------------------------------------
if [ "$0" != './eigen_sparse.sh' ]
then
	echo './eigen_sparse.sh: must be executed from its directory'
	exit 1
fi
# ---------------------------------------------------------------------------
if [ ! -e build ]
then
	mkdir build
fi
cd build
# ---------------------------------------------------------------------------
cat << EOF > eigen_sparse.cpp
# include <iostream>
# include <Eigen/Sparse>
int main(void)
{   typedef Eigen::SparseMatrix<double, Eigen::ColMajor>       sparse_matrix;
	typedef Eigen::SimplicialLDLT<sparse_matrix, Eigen::Lower> sparse_cholesky;
	typedef sparse_matrix::InnerIterator                       column_itr;
	//
	int n = 5;
	sparse_matrix H(n, n), x(n, 1), y(n, 1);
	for(int i = 0; i < n; i++)
	{	H.insert(i, i) = 1.0;         // n by n identity matrix
		x.insert(i, 0) = double(i);   // transpose of (0, 1, ..., n-1)
	}
	y = H * x;
	std::cout << "Test if Eignes's sparse matrix multiply stores zeros\\n";
	for(int j = 0; j < y.outerSize(); j++)
	{	for(column_itr itr(y, j); itr; ++itr)
		{	int row      = itr.row();
			int col      = itr.col();
			double value = itr.value();
			std::cout << "y(" << row << "," << col << ")=" << value << "\\n";
		}
	}
	sparse_cholesky C;  // begin sparse Cholesky solve x = H * y
	C.analyzePattern(H);
	C.factorize(H);
	y = C.solve(x);     // end sparse Cholesky solve
	std::cout << "Test if Eignes's sparse Cholesky solve stores zeros\\n";
	for(int j = 0; j < y.outerSize(); j++)
	{	for(column_itr itr(y, j); itr; ++itr)
		{	int row      = itr.row();
			int col      = itr.col();
			double value = itr.value();
			std::cout << "y(" << row << "," << col << ")=" << value << "\\n";
		}
	}
	return 0;
}
EOF
# ---------------------------------------------------------------------------
echo_eval g++ -I$HOME/prefix/eigen/include eigen_sparse.cpp -o eigen_sparse
echo_eval ./eigen_sparse
