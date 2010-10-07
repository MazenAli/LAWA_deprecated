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


template <typename T, typename Index>
T
NoPreconditioner<T, Index>::operator()(const Index &index) const
{
    return 1.;
}



template <typename T, typename Basis, typename BilinearForm>
DiagonalMatrixPreconditioner1D<T,Basis,BilinearForm>::DiagonalMatrixPreconditioner1D(const BilinearForm &_a)
    : a(_a)
{
}

template <typename T, typename Basis, typename BilinearForm>
T
DiagonalMatrixPreconditioner1D<T,Basis,BilinearForm>::operator()(XType xtype, int j, int k) const
{
    return 1./std::sqrt(fabs(a(xtype,j,k,xtype,j,k)));
}

template <typename T, typename Basis, typename BilinearForm>
T
DiagonalMatrixPreconditioner1D<T,Basis,BilinearForm>::operator()(const Index1D &index) const
{
    return DiagonalMatrixPreconditioner1D<T,Basis,BilinearForm>::operator()(index.xtype, index.j, index.k);
}


template <typename T, typename Basis, typename BilinearForm>
H1Preconditioner1D<T,Basis,BilinearForm>::H1Preconditioner1D(const BilinearForm &_a)
    : basis(_a.getBasis()), phi(basis.mra), psi(basis), d_phi(basis.mra, 1), d_psi(basis, 1),
      integral_sfsf(phi, phi), dd_integral_sfsf(d_phi, d_phi),
      integral_ww(psi,psi), dd_integral_ww(d_psi,d_psi)
{
}

template <typename T, typename Basis, typename BilinearForm>
T
H1Preconditioner1D<T,Basis,BilinearForm>::operator()(XType xtype, int j, int k) const
{
    T val, dd_val;
    if(xtype == XBSpline){
        val = integral_sfsf(j, k, j, k);
        dd_val = dd_integral_sfsf(j, k, j, k);
    }
    else {
        val = integral_ww(j, k, j, k);
        dd_val = dd_integral_ww(j, k, j, k);
    }
    return 1./std::sqrt(val+dd_val);
}

template <typename T, typename Basis, typename BilinearForm>
T
H1Preconditioner1D<T,Basis,BilinearForm>::operator()(const Index1D &index) const
{
    return H1Preconditioner1D<T,Basis,BilinearForm>::operator()(index.xtype,index.j,index.k);
}


template <typename T, typename Basis2D, typename BilinearForm>
DiagonalMatrixPreconditioner2D<T,Basis2D,BilinearForm>::DiagonalMatrixPreconditioner2D(const BilinearForm &_a)
    : a(_a)
{
}

template <typename T, typename Basis2D, typename BilinearForm>
T
DiagonalMatrixPreconditioner2D<T,Basis2D,BilinearForm>::operator()(XType xtype1, int j1, int k1,
                                                                   XType xtype2, int j2, int k2) const
{
    return 1./std::sqrt(fabs(a(xtype1,j1,k1, xtype2,j2,k2,  xtype1,j1,k1, xtype2,j2,k2)));
}

template <typename T, typename Basis2D, typename BilinearForm>
T
DiagonalMatrixPreconditioner2D<T,Basis2D,BilinearForm>::operator()(const Index2D &index) const
{
    return DiagonalMatrixPreconditioner2D<T,Basis2D,BilinearForm>::operator()(index.index1.xtype, index.index1.j, index.index1.k,
                                                                              index.index2.xtype, index.index2.j, index.index2.k);
}



template <typename T, typename Basis3D, typename BilinearForm>
DiagonalMatrixPreconditioner3D<T,Basis3D,BilinearForm>::DiagonalMatrixPreconditioner3D(const BilinearForm &_a)
    : a(_a)
{
}

template <typename T, typename Basis3D, typename BilinearForm>
T
DiagonalMatrixPreconditioner3D<T,Basis3D,BilinearForm>::operator()(XType xtype1, int j1, int k1,
                                                                   XType xtype2, int j2, int k2,
                                                                   XType xtype3, int j3, int k3) const
{
    return 1./std::sqrt(fabs(a(xtype1,j1,k1, xtype2,j2,k2, xtype3,j3,k3,
    	                       xtype1,j1,k1, xtype2,j2,k2, xtype3,j3,k3)));
}

template <typename T, typename Basis3D, typename BilinearForm>
T
DiagonalMatrixPreconditioner3D<T,Basis3D,BilinearForm>::operator()(const Index3D &index) const
{
    return DiagonalMatrixPreconditioner3D<T,Basis3D,BilinearForm>::operator()(index.index1.xtype, index.index1.j, index.index1.k,
                                                                              index.index2.xtype, index.index2.j, index.index2.k,
                                                                              index.index3.xtype, index.index3.j, index.index3.k);
}


template <typename T, typename Basis2D, typename BilinearForm>
SpaceTimePreconditioner2D<T,Basis2D,BilinearForm>::SpaceTimePreconditioner2D(const BilinearForm &_a)
    : a(_a), basis(_a.getBasis()), phi_t(basis.first.mra), phi_x(basis.second.mra),
     psi_t(basis.first), psi_x(basis.second),
     integral_t_sfsf(phi_t, phi_t), integral_x_sfsf(phi_x, phi_x),
     integral_t_ww(psi_t, psi_t), integral_x_ww(psi_x, psi_x)
{
}

template <typename T, typename Basis2D, typename BilinearForm>
T
SpaceTimePreconditioner2D<T,Basis2D,BilinearForm>::operator()(XType xtype1, int j1, int k1,
                                                              XType xtype2, int j2, int k2) const
{
    T val_id_t, val_id_x;
    if(xtype1 == XBSpline){
        val_id_t = integral_t_sfsf(j1, k1, j1, k1);
        if(xtype2 == XBSpline){
            val_id_x = integral_x_sfsf(j2, k2, j2, k2);
        }
        else{
            val_id_x = integral_x_ww(j2, k2, j2, k2);            
        }
    }
    else{
        val_id_t = integral_t_ww(j1, k1, j1, k1);
        if(xtype2 == XBSpline){
            val_id_x = integral_x_sfsf(j2, k2, j2, k2);
        }
        else{
            val_id_x = integral_x_ww(j2, k2, j2, k2);            
        }
    }
    return 1./std::sqrt(val_id_t * pow2i<T>(2*j2) * val_id_x * fabs(a(xtype1,j1,k1, xtype2,j2,k2,  xtype1,j1,k1, xtype2,j2,k2)));
}

template <typename T, typename Basis2D, typename BilinearForm>
T
SpaceTimePreconditioner2D<T,Basis2D,BilinearForm>::operator()(const Index2D &index) const
{
    return SpaceTimePreconditioner2D<T,Basis2D,BilinearForm>::operator()(index.index1.xtype, index.index1.j, index.index1.k,
                                                                         index.index2.xtype, index.index2.j, index.index2.k);
}

} // namespace lawa
