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

template <typename T, typename Basis>
RHSWithPeaks1D<T,Basis>::RHSWithPeaks1D(const Basis &_basis, T (*_f)(T),
								    const DenseVector<Array<T> > &_singularPoints,
								    const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &_deltas,
								    int order)
	: basis(_basis),  phi(basis.mra), psi(basis),
	  f(_f,_singularPoints), deltas(_deltas),
	  integral_sff(phi, f), integral_wf(psi, f)
{
	integral_sff.quadrature.setOrder(order);
	integral_wf.quadrature.setOrder(order);
}

template <typename T, typename Basis>
T
RHSWithPeaks1D<T,Basis>::operator()(XType xtype, int j, int k) const
{
	T ret = 0.;
	if (xtype == XBSpline)  {
		ret = integral_sff(j,k);
		for (int i=1; i<=deltas.numRows(); ++i) {
			ret += deltas(i,2) * phi(deltas(i,1),j,k);
		}
	}
	else {
		ret = integral_wf(j,k);
		for (int i=1; i<=deltas.numRows(); ++i) {
			ret += deltas(i,2) * psi(deltas(i,1),j,k);
		}
	}

	return ret;
}

template <typename T, typename Basis>
T
RHSWithPeaks1D<T,Basis>::operator()(const Index1D &lambda) const
{
	return RHSWithPeaks1D<T,Basis>::operator()(lambda.xtype, lambda.j, lambda.k);
}

}  //namespace lawa
