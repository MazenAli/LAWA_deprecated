/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009 Sebastian Kestler, Kristina Steih, Mario Rometsch, Alexander Stippler.

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

template <typename T, typename Basis2D>
PDEOperator2D<T, Basis2D>::PDEOperator2D(const Basis2D &_basis, T _c)
    : basis(_basis), diffusion_x(basis.first), reaction_x(basis.first),
      diffusion_y(basis.second), reaction_y(basis.second), c(_c)
{
}

template <typename T, typename Basis2D>
T
PDEOperator2D<T,Basis2D>::getc() const
{
    return c;
}

template <typename T, typename Basis2D>
const Basis2D&
PDEOperator2D<T,Basis2D>::getBasis() const
{
    return basis;
}

template <typename T, typename Basis2D>
T
PDEOperator2D<T, Basis2D>::operator()(XType row_xtype_x, int j1_x, int k1_x,
                                      XType row_xtype_y, int j1_y, int k1_y,
                                      XType col_xtype_x, int j2_x, int k2_x,
                                      XType col_xtype_y, int j2_y, int k2_y) const
{
	return    diffusion_x(row_xtype_x,j1_x,k1_x, col_xtype_x, j2_x, k2_x) *
		       reaction_y(row_xtype_y,j1_y,k1_y, col_xtype_y, j2_y, k2_y)
		   +   reaction_x(row_xtype_x,j1_x,k1_x, col_xtype_x, j2_x, k2_x) *
		      diffusion_y(row_xtype_y,j1_y,k1_y, col_xtype_y, j2_y, k2_y)
		   + c*reaction_x(row_xtype_x,j1_x,k1_x, col_xtype_x, j2_x, k2_x) *
		       reaction_y(row_xtype_y,j1_y,k1_y, col_xtype_y, j2_y, k2_y);
}

template <typename T, typename Basis2D>
T
PDEOperator2D<T, Basis2D>::operator()(const Index2D &row_index, const Index2D &col_index) const
{
    return PDEOperator2D<T, Basis2D>::operator()(row_index.index1.xtype, row_index.index1.j, row_index.index1.k,
                                                       row_index.index2.xtype, row_index.index2.j, row_index.index2.k,
                                                       col_index.index1.xtype, col_index.index1.j, col_index.index1.k,
                                                       col_index.index2.xtype, col_index.index2.j, col_index.index2.k);
}

}    //namespace lawa
