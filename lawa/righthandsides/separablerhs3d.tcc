/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Kristina Steih, Alexander Stippler.

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

template<typename T, typename Basis3D>
SeparableRHS3D<T, Basis3D>::SeparableRHS3D(const Basis3D& _basis, const SeparableFunction3D<T>& _F,
		const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &_deltas_x,
		const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &_deltas_y,
		const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &_deltas_z, int order)
    : basis(_basis), F(_F),
      deltas_x(_deltas_x), deltas_y(_deltas_y), deltas_z(_deltas_z),
	  phi_x(_basis.first.mra), phi_y(_basis.second.mra), phi_z(_basis.third.mra),
	  psi_x(_basis.first), psi_y(_basis.second), psi_z(_basis.third),
      integral_sff_x(phi_x, _F.F_x), integral_sff_y(phi_y, _F.F_y), integral_sff_z(phi_z, _F.F_z),
      integral_wf_x(psi_x, _F.F_x),  integral_wf_y(psi_y, _F.F_y), integral_wf_z(psi_z, _F.F_z)
{
	integral_sff_x.quadrature.setOrder(order);
	integral_sff_y.quadrature.setOrder(order);
	integral_sff_z.quadrature.setOrder(order);
	integral_wf_x.quadrature.setOrder(order);
	integral_wf_y.quadrature.setOrder(order);
	integral_wf_z.quadrature.setOrder(order);
}

template<typename T, typename Basis3D>
T
SeparableRHS3D<T, Basis3D>::operator()(XType xtype_x, int j_x, int k_x,
								       XType xtype_y, int j_y, int k_y,
								       XType xtype_z, int j_z, int k_z) const
{
    T val_x = 0;
    T val_y = 0;
    T val_z = 0;

    if(xtype_x == XBSpline){
        val_x = integral_sff_x(j_x, k_x);
        for (int i=1; i<=deltas_x.numRows(); ++i) {
        	val_x += deltas_x(i,2) * phi_x(deltas_x(i,1),j_x,k_x);
        }
    }
    else{
        val_x = integral_wf_x(j_x, k_x);
        for (int i=1; i<=deltas_x.numRows(); ++i) {
        	val_x += deltas_x(i,2) * psi_x(deltas_x(i,1),j_x,k_x);
        }
    }

    if(xtype_y == XBSpline){
        val_y = integral_sff_y(j_y, k_y);
        for (int i=1; i<=deltas_y.numRows(); ++i) {
            val_y += deltas_y(i,2) * phi_y(deltas_y(i,1),j_y,k_y);
        }
    }
    else{
        val_y = integral_wf_y(j_y, k_y);
        for (int i=1; i<=deltas_y.numRows(); ++i) {
             val_y += deltas_y(i,2) * psi_y(deltas_y(i,1),j_y,k_y);
        }
    }

    if(xtype_z == XBSpline){
        val_z = integral_sff_z(j_z, k_z);
        for (int i=1; i<=deltas_z.numRows(); ++i) {
            val_z += deltas_z(i,2) * phi_z(deltas_z(i,1),j_z,k_z);
        }
    }
    else{
        val_z = integral_wf_z(j_z, k_z);
        for (int i=1; i<=deltas_z.numRows(); ++i) {
            val_z += deltas_z(i,2) * psi_z(deltas_z(i,1),j_z,k_z);
        }
    }

    return val_x * val_y * val_z;
}

template<typename T, typename Basis3D>
T
SeparableRHS3D<T, Basis3D>::operator()(const Index3D &index) const
{
    return SeparableRHS3D<T, Basis3D>::operator()(index.index1.xtype, index.index1.j, index.index1.k,
                                                  index.index2.xtype, index.index2.j, index.index2.k,
                                                  index.index3.xtype, index.index3.j, index.index3.k);
}

} // namespace lawa

