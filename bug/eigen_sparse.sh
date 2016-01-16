# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-16 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
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
{	using std::cout;
	using std::endl;
	//
	int n = 5;
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor> sparse_matrix;
	typedef sparse_matrix::InnerIterator col_itr;
	//
	sparse_matrix A(n, n), x(n, 1), y(n, 1);
	for(int i = 0; i < n; i++)
	{	A.insert(i, i) = 1.0;
		x.insert(i, 0) = double(i);
	}
	y = A * x;
	//
	cout << "Test if Eignes's sparse matrices store elements with value zero";
	cout << endl;
	for(int j = 0; j < y.outerSize(); j++)
	{	for(col_itr itr(y, j); itr; ++itr)
		{	int row      = itr.row();
			int col      = itr.col();
			double value = itr.value();
			cout << "y(" << row << "," << col << ")=" << value << endl;
		}
	}
	return 0;
}
EOF
# ---------------------------------------------------------------------------
echo_eval g++ -I$HOME/prefix/eigen/include eigen_sparse.cpp -o eigen_sparse
echo_eval ./eigen_sparse

