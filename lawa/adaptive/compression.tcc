/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009 Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

namespace lawa {

template <typename T, typename Index, typename Basis>
NoCompression<T,Index,Basis>::NoCompression(const Basis &_basis)
    : basis(_basis)
{
}

template <typename T, typename Index, typename Basis>
void
NoCompression<T,Index,Basis>::setParameters(const IndexSet<Index> &LambdaRow)
{

}

template <typename T, typename Index, typename Basis>
IndexSet<Index>
NoCompression<T,Index,Basis>::SparsityPattern(const Index &lambda_col, const IndexSet<Index> &LambdaRow)
{
	return LambdaRow;
}

template <typename T, typename Basis>
CompressionPDE1D<T,Basis>::CompressionPDE1D(const Basis &_basis)
    : basis(_basis), s_tilde(-1), jmin(100), jmax(-30)
{
}

template <typename T, typename Basis>
void
CompressionPDE1D<T,Basis>::setParameters(const IndexSet<Index1D> &LambdaRow) {
	typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;
	s_tilde = -1;
	jmin = 100;
	jmax = -30;
	for (set1d_const_it lambda_row = LambdaRow.begin(); lambda_row != LambdaRow.end(); ++lambda_row) {
	    jmin = std::min(jmin,(*lambda_row).j);
	    jmax = std::max(jmax,(*lambda_row).j);
	}
	s_tilde = jmax-jmin;
}

template <typename T, typename Basis>
IndexSet<Index1D>
CompressionPDE1D<T,Basis>::SparsityPattern(const Index1D &lambda_col, const IndexSet<Index1D> &LambdaRow) {
	IndexSet<Index1D> LambdaRowSparse(LambdaRow.d,LambdaRow.d_);
	LambdaRowSparse = (lambda_col, basis, s_tilde, jmin, jmax, false);
	return LambdaRowSparse;
}


template <typename T, typename Basis>
CompressionPDE2D<T,Basis>::CompressionPDE2D(const Basis &_basis)
    : basis(_basis), s_tilde_x(-1), jmin_x(100), jmax_x(-30), s_tilde_y(-1), jmin_y(100), jmax_y(-30)
{
}

template <typename T, typename Basis>
void
CompressionPDE2D<T,Basis>::setParameters(const IndexSet<Index2D> &LambdaRow) {
	typedef typename IndexSet<Index2D>::const_iterator set2d_const_it;
	jmin_x = 100, jmax_x=-30, jmin_y = 100, jmax_y=-30;
	for (set2d_const_it lambda_col = LambdaRow.begin(); lambda_col != LambdaRow.end(); ++lambda_col) {
	    jmin_x = std::min(jmin_x,(*lambda_col).index1.j);
	    jmax_x = std::max(jmax_x,(*lambda_col).index1.j);
	    jmin_y = std::min(jmin_y,(*lambda_col).index2.j);
	    jmax_y = std::max(jmax_y,(*lambda_col).index2.j);
	}
	s_tilde_x = jmax_x-jmin_x;
	s_tilde_y = jmax_y-jmin_y;
}

template <typename T, typename Basis>
IndexSet<Index2D>
CompressionPDE2D<T,Basis>::SparsityPattern(const Index2D &lambda_col, const IndexSet<Index2D> &LambdaRow) {
	typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;
	IndexSet<Index2D> LambdaRowSparse(LambdaRow.d,LambdaRow.d_);

	IndexSet<Index1D> Lambda_x = lambdaTilde1d_PDE(lambda_col.index1, basis.first, s_tilde_x, jmin_x, jmax_x, false);
	IndexSet<Index1D> Lambda_y = lambdaTilde1d_PDE(lambda_col.index2, basis.second, s_tilde_y, jmin_y, jmax_y, false);
	for (set1d_const_it lambda_x = Lambda_x.begin(); lambda_x != Lambda_x.end(); ++lambda_x) {
		for (set1d_const_it lambda_y = Lambda_y.begin(); lambda_y != Lambda_y.end(); ++lambda_y) {
			Index2D index2d(*lambda_x,*lambda_y);
			if (LambdaRow.count(index2d) > 0) {
				LambdaRowSparse.insert(index2d);
			}
		}
	}
	return LambdaRowSparse;
}

template <typename T, typename Basis>
IndexSet<Index2D>
CompressionPDE2D<T,Basis>::SparsityPattern(const Index2D &lambda_col, int jmin_x, int jmin_y, int s_tilde) {
	typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;
	IndexSet<Index2D> ret(basis.first.d,basis.first.d_);

	T factor;
	if      (basis.first.d == 2) factor = 2.;
	else if (basis.first.d == 3) factor = 2./3.;
	else {
		std::cout << "CompressionPDE2D not implemented for orders higher than 3." << std::endl;
		exit(1);
	}

	for (int s_tilde_x=0; s_tilde_x <= factor*s_tilde; ++s_tilde_x) {
		for (int s_tilde_y=0; s_tilde_y <= factor*s_tilde; ++s_tilde_y) {
			if (s_tilde_x+s_tilde_y <= factor*s_tilde) {
				IndexSet<Index1D> Lambda_x = lambdaTilde1d_PDE(lambda_col.index1, basis.first,  s_tilde_x, jmin_x,
														       lambda_col.index1.j+s_tilde_x, false);
			    IndexSet<Index1D> Lambda_y = lambdaTilde1d_PDE(lambda_col.index2, basis.second, s_tilde_y, jmin_y,
					                                            lambda_col.index2.j+s_tilde_y, false);


			    for (set1d_const_it lambda_x = Lambda_x.begin(); lambda_x != Lambda_x.end(); ++lambda_x) {
				    for (set1d_const_it lambda_y = Lambda_y.begin(); lambda_y != Lambda_y.end(); ++lambda_y) {
					    ret.insert(Index2D(*lambda_x,*lambda_y));
				    }
			    }
			}
		}
	}

	return ret;
}


}	//namespace lawa
