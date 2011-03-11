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

template <typename T, typename Basis3D>
HelmholtzOperator3D<T, Basis3D>::HelmholtzOperator3D(const Basis3D & _basis, const T _c)
    : basis(_basis), c(_c),
      dd_x(basis.first), id_x(basis.first), dd_y(basis.second), id_y(basis.second),
      dd_z(basis.third), id_z(basis.third),
      phi_x(basis.first.mra), d_phi_x(basis.first.mra, 1),
      phi_y(basis.second.mra), d_phi_y(basis.second.mra, 1),
      phi_z(basis.third.mra), d_phi_z(basis.third.mra, 1),
      psi_x(basis.first), d_psi_x(basis.first, 1),
      psi_y(basis.second), d_psi_y(basis.second, 1),
      psi_z(basis.third), d_psi_z(basis.third, 1),
      integral_sfsf_x(phi_x, phi_x), dd_integral_sfsf_x(d_phi_x, d_phi_x),
      integral_sfsf_y(phi_y, phi_y), dd_integral_sfsf_y(d_phi_y, d_phi_y),
      integral_sfsf_z(phi_z, phi_z), dd_integral_sfsf_z(d_phi_z, d_phi_z),
      integral_sfw_x(phi_x, psi_x), dd_integral_sfw_x(d_phi_x, d_psi_x),
      integral_sfw_y(phi_y, psi_y), dd_integral_sfw_y(d_phi_y, d_psi_y),
      integral_sfw_z(phi_z, psi_z), dd_integral_sfw_z(d_phi_z, d_psi_z),
      integral_wsf_x(psi_x, phi_x), dd_integral_wsf_x(d_psi_x, d_phi_x),
      integral_wsf_y(psi_y, phi_y), dd_integral_wsf_y(d_psi_y, d_phi_y),
      integral_wsf_z(psi_z, phi_z), dd_integral_wsf_z(d_psi_z, d_phi_z),
      integral_ww_x(psi_x, psi_x), dd_integral_ww_x(d_psi_x, d_psi_x),
      integral_ww_y(psi_y, psi_y), dd_integral_ww_y(d_psi_y, d_psi_y),
      integral_ww_z(psi_z, psi_z), dd_integral_ww_z(d_psi_z, d_psi_z)
{
}

template <typename T, typename Basis3D>
T
HelmholtzOperator3D<T,Basis3D>::getc() const
{
    return c;
}

template <typename T, typename Basis3D>
const Basis3D&
HelmholtzOperator3D<T,Basis3D>::getBasis() const
{
    return basis;
}

template <typename T, typename Basis3D>
T
HelmholtzOperator3D<T, Basis3D>::operator()(XType row_xtype_x, int j1_x, int k1_x,
                                            XType row_xtype_y, int j1_y, int k1_y,
                                            XType row_xtype_z, int j1_z, int k1_z,
                                            XType col_xtype_x, int j2_x, int k2_x,
                                            XType col_xtype_y, int j2_y, int k2_y,
                                            XType col_xtype_z, int j2_z, int k2_z) const
{
    T val_x = 0;
    T dd_val_x = 0;
    T val_y = 0;
    T dd_val_y = 0;
    T val_z = 0;
    T dd_val_z = 0;

    if(row_xtype_x == XBSpline){
         if(col_xtype_x == XBSpline){
             val_x =       integral_sfsf_x(j1_x, k1_x, j2_x, k2_x);
             dd_val_x = dd_integral_sfsf_x(j1_x, k1_x, j2_x, k2_x);
         }
         else{
             val_x =       integral_sfw_x(j1_x, k1_x, j2_x, k2_x);
             dd_val_x = dd_integral_sfw_x(j1_x, k1_x, j2_x, k2_x);
         }
    }
    else{
        if(col_xtype_x == XBSpline){
            val_x = integral_wsf_x(j1_x, k1_x, j2_x, k2_x);
            dd_val_x = dd_integral_wsf_x(j1_x, k1_x, j2_x, k2_x);
        }
        else{
            val_x =       integral_ww_x(j1_x, k1_x, j2_x, k2_x);
            dd_val_x = dd_integral_ww_x(j1_x, k1_x, j2_x, k2_x);
        }
    }

    if(row_xtype_y == XBSpline){
         if(col_xtype_y == XBSpline){
             val_y =       integral_sfsf_y(j1_y, k1_y, j2_y, k2_y);
             dd_val_y = dd_integral_sfsf_y(j1_y, k1_y, j2_y, k2_y);
         }
         else{
             val_y = integral_sfw_y(j1_y, k1_y, j2_y, k2_y);
             dd_val_y = dd_integral_sfw_y(j1_y, k1_y, j2_y, k2_y);
         }
    }
    else{
         if(col_xtype_y == XBSpline){
             val_y = integral_wsf_y(j1_y, k1_y, j2_y, k2_y);
             dd_val_y = dd_integral_wsf_y(j1_y, k1_y, j2_y, k2_y);
         }
         else{
             val_y = integral_ww_y(j1_y, k1_y, j2_y, k2_y);
             dd_val_y = dd_integral_ww_y(j1_y, k1_y, j2_y, k2_y);
         }
    }

    if(row_xtype_z == XBSpline){
        if(col_xtype_z == XBSpline){
            val_z =       integral_sfsf_z(j1_z, k1_z, j2_z, k2_z);
            dd_val_z = dd_integral_sfsf_z(j1_z, k1_z, j2_z, k2_z);
        }
        else{
            val_z = integral_sfw_z(j1_z, k1_z, j2_z, k2_z);
            dd_val_z = dd_integral_sfw_z(j1_z, k1_z, j2_z, k2_z);
        }
    }
    else{
        if(col_xtype_z == XBSpline){
            val_z = integral_wsf_z(j1_z, k1_z, j2_z, k2_z);
            dd_val_z = dd_integral_wsf_z(j1_z, k1_z, j2_z, k2_z);
        }
        else{
            val_z = integral_ww_z(j1_z, k1_z, j2_z, k2_z);
            dd_val_z = dd_integral_ww_z(j1_z, k1_z, j2_z, k2_z);
        }
    }

    return dd_val_x * val_y * val_z + val_x * dd_val_y * val_z + val_x * val_y * dd_val_z + c * val_x * val_y * val_z;
}

template <typename T, typename Basis3D>
T
HelmholtzOperator3D<T, Basis3D>::operator()(const Index3D &row_index, const Index3D &col_index) const
{
    return HelmholtzOperator3D<T, Basis3D>::operator()(row_index.index1.xtype, row_index.index1.j, row_index.index1.k,
                                                       row_index.index2.xtype, row_index.index2.j, row_index.index2.k,
                                                       row_index.index3.xtype, row_index.index3.j, row_index.index3.k,
                                                       col_index.index1.xtype, col_index.index1.j, col_index.index1.k,
                                                       col_index.index2.xtype, col_index.index2.j, col_index.index2.k,
                                                       col_index.index3.xtype, col_index.index3.j, col_index.index3.k);
}

}    //namespace lawa
