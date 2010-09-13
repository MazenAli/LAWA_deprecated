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

template<typename T, typename Basis2D, QuadratureType Quad>
SmoothRHSWithAlignedSing2D<T,Basis2D,Quad>::SmoothRHSWithAlignedSing2D(const Basis2D &_basis,
												const Function2D<T> &_F, int order)
   : basis(_basis),
     phi_x(_basis.first.mra), phi_y(_basis.second.mra),
     psi_x(_basis.first), psi_y(_basis.second),
     integral2d_phi_x_phi_y(_F, phi_x, phi_y),
     integral2d_psi_x_phi_y(_F, psi_x, phi_y),
     integral2d_phi_x_psi_y(_F, phi_x, psi_y),
     integral2d_psi_x_psi_y(_F, psi_x, psi_y)
{
	integral2d_phi_x_phi_y.quadrature.setOrder(order);
	integral2d_psi_x_phi_y.quadrature.setOrder(order);
	integral2d_phi_x_psi_y.quadrature.setOrder(order);
	integral2d_psi_x_psi_y.quadrature.setOrder(order);
}

template<typename T, typename Basis2D, QuadratureType Quad>
SmoothRHSWithAlignedSing2D<T,Basis2D,Quad>::SmoothRHSWithAlignedSing2D(const Basis2D &_basis,
												const Function2D<T> &_F, int order,
												int deriv_x, int deriv_y)
   : basis(_basis),
     phi_x(_basis.first.mra, deriv_x), phi_y(_basis.second.mra, deriv_y),
     psi_x(_basis.first, deriv_x), psi_y(_basis.second, deriv_y),
     integral2d_phi_x_phi_y(_F, phi_x, phi_y),
     integral2d_psi_x_phi_y(_F, psi_x, phi_y),
     integral2d_phi_x_psi_y(_F, phi_x, psi_y),
     integral2d_psi_x_psi_y(_F, psi_x, psi_y)
{
	integral2d_phi_x_phi_y.quadrature.setOrder(order);
	integral2d_psi_x_phi_y.quadrature.setOrder(order);
	integral2d_phi_x_psi_y.quadrature.setOrder(order);
	integral2d_psi_x_psi_y.quadrature.setOrder(order);
}

template<typename T, typename Basis2D, QuadratureType Quad>
T
SmoothRHSWithAlignedSing2D<T,Basis2D,Quad>::operator()(XType xtype_x, int j_x, int k_x,
                                                     XType xtype_y, int j_y, int k_y) const
{
	if(xtype_x == XBSpline){
		if(xtype_y == XBSpline){
			return integral2d_phi_x_phi_y(j_x,k_x,j_y,k_y);
		}
		else {
			return integral2d_phi_x_psi_y(j_x,k_x,j_y,k_y);
		}
	}
	else {
		if(xtype_y == XBSpline){
			return integral2d_psi_x_phi_y(j_x,k_x,j_y,k_y);
		}
		else {
			return integral2d_psi_x_psi_y(j_x,k_x,j_y,k_y);
		}
	}
}

template<typename T, typename Basis2D, QuadratureType Quad>
T
SmoothRHSWithAlignedSing2D<T,Basis2D,Quad>::operator()(const Index2D &index) const
{
	return SmoothRHSWithAlignedSing2D<T,Basis2D,Quad>::operator()(index.index1.xtype, index.index1.j, index.index1.k,
	                                                            index.index2.xtype, index.index2.j, index.index2.k);
}

}	//namespace lawa
