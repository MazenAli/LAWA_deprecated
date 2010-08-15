/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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

template <typename T, typename Basis, typename Preconditioner>
RHS<T,Index1d,Basis,Preconditioner>::RHS(const Basis &_basis,T (*_f)(T), const DenseVector<Array<T> > &_singularPoints, GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > _deltas)
	: basis(_basis),  phi(basis.mra), psi(basis),
	  f(_f,_singularPoints), deltas(_deltas),
	  integral_sff(phi, f), integral_wf(psi, f),
	  P(), rhs(basis.d_,basis.d_), rhs_abs(basis.d_,basis.d_)
{

}

template <typename T, typename Basis, typename Preconditioner>
Coefficients<Lexicographical,T,Index1d>
RHS<T,Index1d,Basis,Preconditioner>::operator()(const IndexSet<Index1d> &Lambda) {
	typedef typename IndexSet<Index1d>::iterator const_set_it;
	Coefficients<Lexicographical,T,Index1d> ret(Lambda.d,Lambda.d_);
	for (const_set_it lambda = Lambda.begin(); lambda != Lambda.end(); ++lambda) {
		if (rhs.count((*lambda)) == 0) {
			T tmp = 0.;
			int j=(*lambda).j, k=(*lambda).k;
			if ((*lambda).xtype == XBSpline)  {
				tmp = integral_sff(j,k);
				for (int i=1; i<=deltas.numRows(); ++i) tmp += deltas(i,1) * phi(deltas(i,2),j,k);
				tmp *= P(*lambda);
			}
			else {
				tmp = integral_wf(j,k);
				for (int i=1; i<=deltas.numRows(); ++i) tmp += deltas(i,1) * psi(deltas(i,2),j,k);
				tmp *= P(*lambda);
			}
			ret[*lambda] = tmp;
			rhs[*lambda] = tmp;
		}
		else {
			ret[*lambda] = rhs[*lambda];
		}
	}
	return ret;
}

}  //namespace lawa
