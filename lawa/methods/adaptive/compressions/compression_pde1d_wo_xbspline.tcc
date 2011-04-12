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

template <typename T, typename Index, typename Basis1D>
CompressionPDE1D_WO_XBSpline<T,Index,Basis1D>::CompressionPDE1D_WO_XBSpline(const Basis1D &_basis)
    : basis(_basis), s_tilde(0), jmin(100), jmax(-30)
{
}

template <typename T, typename Index, typename Basis1D>
void
CompressionPDE1D_WO_XBSpline<T,Index,Basis1D>::setParameters(const IndexSet<Index> &Lambda) {
	typedef typename IndexSet<Index>::const_iterator set1d_const_it;
	s_tilde = -1;
	jmin = 100;
	jmax = -30;
	for (set1d_const_it lambda = Lambda.begin(); lambda != Lambda.end(); ++lambda) {
	    jmin = std::min(jmin,(*lambda).j);
	    jmax = std::max(jmax,(*lambda).j);
	}
	s_tilde = jmax-jmin;
}

template <typename T, typename Index, typename Basis1D>
IndexSet<Index>
CompressionPDE1D_WO_XBSpline<T,Index,Basis1D>::SparsityPattern(const Index &lambda,
															   const IndexSet<Index> &Lambda, int J) {
	typedef typename IndexSet<Index>::const_iterator set1d_const_it;

	IndexSet<Index> LambdaRowSparse(Lambda.d,Lambda.d_), Lambda_x(Lambda.d,Lambda.d_);
	if (J==-1) {
		Lambda_x = lambdaTilde1d_PDE_WO_XBSpline(lambda, basis, s_tilde, jmin, jmax);
	}
	else {
		Lambda_x = lambdaTilde1d_PDE_WO_XBSpline(lambda, basis, std::min(int(s_tilde),J), jmin, jmax);
	}

	for (set1d_const_it lambda_x = Lambda_x.begin(); lambda_x != Lambda_x.end(); ++lambda_x) {
		if (Lambda.count(*lambda_x)>0) {
			LambdaRowSparse.insert(*lambda_x);
		}
	}
	return LambdaRowSparse;
}

}	//namespace lawa
