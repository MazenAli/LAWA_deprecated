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


template <typename T>
T
NoPreconditioner1D<T>::operator()(XType xtype, int j, int k) const
{
    return 1.;
}

template <typename T>
T
NoPreconditioner1D<T>::operator()(const Index1D &index) const
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

} // namespace lawa
