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
	typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;

	IndexSet<Index1D> LambdaRowSparse(LambdaRow.d,LambdaRow.d_);
	int s = std::max(abs(lambda_col.j-jmin),abs(lambda_col.j-jmax));
	s = std::max(s,int(s_tilde));
	IndexSet<Index1D> Lambda_x = lambdaTilde1d_PDE(lambda_col, basis, s, jmin, jmax, false);
	for (set1d_const_it lambda_x = Lambda_x.begin(); lambda_x != Lambda_x.end(); ++lambda_x) {
		if (LambdaRow.count(*lambda_x)>0) {
			LambdaRowSparse.insert(*lambda_x);
		}
	}
	return LambdaRowSparse;
}


template <typename T>
CompressionCGMYOperator1D<T,Basis<T,Primal,R,CDF> >::CompressionCGMYOperator1D(const Basis<T,Primal,R,CDF> &_basis, T Y)
    : basis(_basis), jmin(100), jmax(-30), compr_c(1.), psi(basis)
{
	compr_c = T(2*basis.d)/T(2*basis.d_+Y);
}

template <typename T>
void
CompressionCGMYOperator1D<T,Basis<T,Primal,R,CDF> >::setParameters(const IndexSet<Index1D> &LambdaRow) {
	typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;
	jmin = 100;
	jmax = -30;
	for (set1d_const_it lambda_row = LambdaRow.begin(); lambda_row != LambdaRow.end(); ++lambda_row) {
	    jmin = std::min(jmin,(*lambda_row).j);
	    jmax = std::max(jmax,(*lambda_row).j);
	}
}

template <typename T>
IndexSet<Index1D>
CompressionCGMYOperator1D<T,Basis<T,Primal,R,CDF> >::SparsityPattern(const Index1D &lambda_col, const IndexSet<Index1D> &LambdaRow) {
	typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;

	IndexSet<Index1D> LambdaRowSparse(LambdaRow.d,LambdaRow.d_);
	if (lambda_col.xtype == XBSpline) {
		return LambdaRow;
	}
	for (set1d_const_it lambda = LambdaRow.begin(); lambda!=LambdaRow.end(); ++lambda) {
		if ((*lambda).xtype == XBSpline) {
			LambdaRowSparse.insert(*lambda);
		}
		else {
			T delta = 3*compr_c*std::max(pow2i<T>(-std::min(lambda_col.j,(*lambda).j)),
									     pow2i<T>(-(jmax-1) + compr_alpha*(2*(jmax-1)-lambda_col.j-(*lambda).j)));
			if (lawa::distance(psi.support(lambda_col.j,lambda_col.k),psi.support((*lambda).j,(*lambda).k)) < delta) {
				LambdaRowSparse.insert(*lambda);
			}
		}
	}
	return LambdaRowSparse;
}


template <typename T, typename Compression_x, typename Compression_y>
TensorCompression2D<T,Compression_x,Compression_y>::TensorCompression2D(Compression_x &_c_x,
																	    Compression_y &_c_y)
    : c_x(_c_x), c_y(_c_y)
{
}

template <typename T, typename Compression_x, typename Compression_y>
void
TensorCompression2D<T,Compression_x,Compression_y>::setParameters(const IndexSet<Index2D> &LambdaRow)
{
	int d=LambdaRow.d, d_=LambdaRow.d_;
	IndexSet<Index1D> LambdaRow_x(d,d_), LambdaRow_y(d,d_);
	split(LambdaRow, LambdaRow_x, LambdaRow_y);
	c_x.setParameters(LambdaRow_x);
	c_y.setParameters(LambdaRow_y);
}

template <typename T, typename Compression_x, typename Compression_y>
IndexSet<Index2D>
TensorCompression2D<T,Compression_x,Compression_y>::SparsityPattern(const Index2D &lambda_col,
																	const IndexSet<Index2D> &LambdaRow)
{
	typedef typename IndexSet<Index2D>::const_iterator set2d_const_it;

	int d=LambdaRow.d, d_=LambdaRow.d_;
	IndexSet<Index1D> LambdaRow_x(d,d_), LambdaRow_y(d,d_);
	split(LambdaRow, LambdaRow_x, LambdaRow_y);
	IndexSet<Index1D> ret_x = c_x.SparsityPattern(lambda_col.index1, LambdaRow_x);
	IndexSet<Index1D> ret_y = c_y.SparsityPattern(lambda_col.index2, LambdaRow_y);

	IndexSet<Index2D> LambdaRowSparse(d,d_);
	for (set2d_const_it lambda=LambdaRow.begin(); lambda!=LambdaRow.end(); ++lambda) {
		if ((LambdaRow_x.count((*lambda).index1)>0) && (LambdaRow_y.count((*lambda).index2)>0))  {
			LambdaRowSparse.insert(*lambda);
		}
	}
	return LambdaRowSparse;
}


template <typename T, typename Basis>
CompressionPDE2D<T,Basis>::CompressionPDE2D(const Basis &_basis, bool _levelthresh)
    : basis(_basis), s_tilde_x(-1), jmin_x(100), jmax_x(-30), s_tilde_y(-1), jmin_y(100), jmax_y(-30), J(0),
      levelthresh(_levelthresh)
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
	if (levelthresh)  {
		J = std::max(s_tilde_x,s_tilde_y);
		J = std::max(10,int(J));
	}
	else			  J = 100;
	std::cout << "   Compression: s_tilde_x = " << s_tilde_x << ", s_tilde_y = " << s_tilde_y << ", J = " << J << std::endl;
}

template <typename T, typename Basis>
IndexSet<Index2D>
CompressionPDE2D<T,Basis>::SparsityPattern(const Index2D &lambda_col, const IndexSet<Index2D> &LambdaRow) {
	typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;
	typedef typename IndexSet<Index2D>::const_iterator set2d_const_it;

	IndexSet<Index2D> LambdaRowSparse(LambdaRow.d,LambdaRow.d_);
	IndexSet<Index1D> Lambda_x = lambdaTilde1d_PDE(lambda_col.index1, basis.first,  s_tilde_x, jmin_x, jmax_x, false);
	IndexSet<Index1D> Lambda_y = lambdaTilde1d_PDE(lambda_col.index2, basis.second, s_tilde_y, jmin_y, jmax_y, false);

	for (set2d_const_it lambda=LambdaRow.begin(); lambda!=LambdaRow.end(); ++lambda) {
		short level_diff = fabs((*lambda).index1.j-lambda_col.index1.j) + fabs((*lambda).index2.j-lambda_col.index2.j);
		if ( (0.5+LambdaRow.d-2)*level_diff > J ) {
			continue;
		}
		if ((Lambda_x.count((*lambda).index1)>0) && (Lambda_y.count((*lambda).index2)>0))  {
			LambdaRowSparse.insert(*lambda);
		}
	}
	return LambdaRowSparse;
}

template <typename T, typename Basis>
IndexSet<Index2D>
CompressionPDE2D<T,Basis>::SparsityPattern(const Index2D &lambda_col, int jmin_x, int jmin_y, int s_tilde,
		                                   int deriv_x, int deriv_y) {
	typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;
	IndexSet<Index2D> ret(basis.first.d,basis.first.d_);

	T factor_x = basis.first.d -  2+1.5-deriv_x;
	T factor_y = basis.second.d - 2+1.5-deriv_y;

	int level_bound_x = round((1./factor_x)*s_tilde);
	int level_bound_y = round((1./factor_y)*s_tilde);

	//std::cout << "factor_x = " << factor_x << ", factor_y = " << factor_y << std::endl;


	for (int s_tilde_x=0; s_tilde_x <= level_bound_x; ++s_tilde_x) {
		for (int s_tilde_y=0; s_tilde_y <= level_bound_y; ++s_tilde_y) {
			if (factor_x*s_tilde_x + factor_y*s_tilde_y <= s_tilde) {
				IndexSet<Index1D> Lambda_x = lambdaTilde1d_PDE(lambda_col.index1, basis.first,  s_tilde_x, jmin_x,
														       lambda_col.index1.j+s_tilde_x, false);
			    IndexSet<Index1D> Lambda_y = lambdaTilde1d_PDE(lambda_col.index2, basis.second, s_tilde_y, jmin_y,
					                                           lambda_col.index2.j+s_tilde_y, false);


			    for (set1d_const_it lambda_x = Lambda_x.begin(); lambda_x != Lambda_x.end(); ++lambda_x) {
				    for (set1d_const_it lambda_y = Lambda_y.begin(); lambda_y != Lambda_y.end(); ++lambda_y) {
				    	if (deriv_x!=0 || deriv_y!=0) {
				    		ret.insert(Index2D(*lambda_x,*lambda_y));
				    	}
				    	else {
				    		int tmp = std::max(lambda_col.index1.j,lambda_col.index2.j)+std::max((*lambda_x).j,(*lambda_y).j);
				    		if (tmp <= s_tilde) {
				    			ret.insert(Index2D(*lambda_x,*lambda_y));
				    		}
				    	}
				    }
			    }
			}
		}
	}

	return ret;
}


template <typename T, typename Basis>
CompressionPDE3D<T,Basis>::CompressionPDE3D(const Basis &_basis)
    : basis(_basis), s_tilde_x(-1), jmin_x(100), jmax_x(-30), s_tilde_y(-1), jmin_y(100), jmax_y(-30),
      s_tilde_z(-1), jmin_z(100), jmax_z(-30)
{
}

template <typename T, typename Basis>
void
CompressionPDE3D<T,Basis>::setParameters(const IndexSet<Index3D> &LambdaRow) {
	typedef typename IndexSet<Index3D>::const_iterator set3d_const_it;
	jmin_x = 100, jmax_x=-30, jmin_y = 100, jmax_y=-30, jmin_z = 100, jmax_z=-30;
	for (set3d_const_it lambda_col = LambdaRow.begin(); lambda_col != LambdaRow.end(); ++lambda_col) {
	    jmin_x = std::min(jmin_x,(*lambda_col).index1.j);
	    jmax_x = std::max(jmax_x,(*lambda_col).index1.j);
	    jmin_y = std::min(jmin_y,(*lambda_col).index2.j);
	    jmax_y = std::max(jmax_y,(*lambda_col).index2.j);
	    jmin_z = std::min(jmin_z,(*lambda_col).index3.j);
	    jmax_z = std::max(jmax_z,(*lambda_col).index3.j);
	}
	s_tilde_x = jmax_x-jmin_x;
	s_tilde_y = jmax_y-jmin_y;
	s_tilde_z = jmax_z-jmin_z;
}

template <typename T, typename Basis>
IndexSet<Index3D>
CompressionPDE3D<T,Basis>::SparsityPattern(const Index3D &lambda_col, const IndexSet<Index3D> &LambdaRow) {
	typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;
	typedef typename IndexSet<Index3D>::const_iterator set3d_const_it;

	set3d_const_it LambdaRow_end = LambdaRow.end();
	IndexSet<Index3D> LambdaRowSparse(LambdaRow.d,LambdaRow.d_);

	IndexSet<Index1D> Lambda_x = lambdaTilde1d_PDE(lambda_col.index1, basis.first, s_tilde_x, jmin_x, jmax_x, false);
	IndexSet<Index1D> Lambda_y = lambdaTilde1d_PDE(lambda_col.index2, basis.second, s_tilde_y, jmin_y, jmax_y, false);
	IndexSet<Index1D> Lambda_z = lambdaTilde1d_PDE(lambda_col.index3, basis.third, s_tilde_z, jmin_z, jmax_z, false);

	for (set3d_const_it it=LambdaRow.begin(); it!=LambdaRow.end(); ++it) {
		if (Lambda_x.count((*it).index1) >0 ) {
			if (Lambda_y.count((*it).index2) >0 ) {
				if (Lambda_z.count((*it).index3) >0 ) {
					LambdaRowSparse.insert(Index3D((*it).index1,(*it).index2,(*it).index3));
				}
			}
		}
	}
/*
	for (set1d_const_it lambda_x = Lambda_x.begin(); lambda_x != Lambda_x.end(); ++lambda_x) {
		for (set1d_const_it lambda_y = Lambda_y.begin(); lambda_y != Lambda_y.end(); ++lambda_y) {
			for (set1d_const_it lambda_z = Lambda_z.begin(); lambda_z != Lambda_z.end(); ++lambda_z) {
				Index3D index3d(*lambda_x,*lambda_y,*lambda_z);
				set3d_const_it it = LambdaRow.find(index3d);
				if (it != LambdaRow_end) {
					LambdaRowSparse.insert(index3d);
				}
				//if (LambdaRow.count(index3d) > 0) {
				//	LambdaRowSparse.insert(index3d);
				//}
			}
		}
	}
*/
	return LambdaRowSparse;
}


}	//namespace lawa
