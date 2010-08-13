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
HelmholtzOperator1d<T,Basis>::HelmholtzOperator1d(const Basis &_basis, const T &_c)
: 	basis(_basis), c(_c),
	phi(basis.mra), d_phi(basis.mra, 1), psi(basis), d_psi(basis, 1),
    integral_sfsf(phi, phi), dd_integral_sfsf(d_phi, d_phi),
    integral_sfw(phi, psi),  dd_integral_sfw(d_phi, d_psi),
    integral_wsf(psi, phi),  dd_integral_wsf(d_psi, d_phi),
    integral_ww(psi, psi),   dd_integral_ww(d_psi, d_psi)
{
}

template <typename T, typename Basis>
HelmholtzOperator1d<T,Basis>::HelmholtzOperator1d(const HelmholtzOperator1d<T,Basis> &_a)
: 	basis(_a.getBasis()), c(_a.getc()),
	phi(basis.mra), d_phi(basis.mra, 1), psi(basis), d_psi(basis, 1),
	integral_sfsf(phi, phi), dd_integral_sfsf(d_phi, d_phi),
	integral_sfw(phi, psi),  dd_integral_sfw(d_phi, d_psi),
	integral_wsf(psi, phi),  dd_integral_wsf(d_psi, d_phi),
	integral_ww(psi, psi),   dd_integral_ww(d_psi, d_psi)
{
}

template <typename T, typename Basis>
T
HelmholtzOperator1d<T,Basis>::getc() const
{
	return c;
}

template <typename T, typename Basis>
const Basis&
HelmholtzOperator1d<T,Basis>::getBasis() const
{
	return basis;
}

template <typename T, typename Basis>
T
HelmholtzOperator1d<T,Basis>::operator()(const Index1d &row_index, const Index1d &col_index) const
{
	T    val = 0.;
	T dd_val = 0.;
	int j_row = row_index.j, k_row = row_index.k, j_col = col_index.j, k_col = col_index.k;

	if ((row_index.xtype == XBSpline) && (col_index.xtype == XBSpline)) {
		dd_val = dd_integral_sfsf(j_row,k_row,j_col,k_col);
		if (c!=0) val = c*integral_sfsf(j_row,k_row,j_col,k_col);
	}
	else if ((row_index.xtype == XBSpline) && (col_index.xtype == XWavelet)) {
		dd_val = dd_integral_sfw(j_row,k_row,j_col,k_col);
		if (c!=0) val = c*integral_sfw(j_row,k_row,j_col,k_col);
	}
	else if ((row_index.xtype == XWavelet) && (col_index.xtype == XBSpline)) {
		dd_val = dd_integral_wsf(j_row,k_row,j_col,k_col);
		if (c!=0) val = c*integral_wsf(j_row,k_row,j_col,k_col);
	}
	else {
		dd_val = dd_integral_ww(j_row,k_row,j_col,k_col);
		if (c!=0) val = c*integral_ww(j_row,k_row,j_col,k_col);
	}
	return val+dd_val;
}

} // namespace lawa
