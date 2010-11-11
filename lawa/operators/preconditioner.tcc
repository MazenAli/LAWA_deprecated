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


template <typename T>
T
DiagonalLevelPreconditioner1D<T>::operator()(XType xtype, int j, int k) const
{
    return 1./std::sqrt(1+pow2i<T>(2*j));
}

template <typename T>
T
DiagonalLevelPreconditioner1D<T>::operator()(const Index1D &index) const
{
	return 1./std::sqrt(1+pow2i<T>(2*index.j));
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
    : basis(_a.getBasis()), phi(basis.mra), d_phi(basis.mra, 1), psi(basis), d_psi(basis, 1),
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


template <typename T, typename Basis2D>
RightNormPreconditioner2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::RightNormPreconditioner2D
	(const SpaceTimeHeatOperator1D<T,Basis2D> &_a)
    : a(_a),
    phi_t(a.basis.first.mra), d_phi_t(a.basis.first.mra, 1),
    phi_x(a.basis.second.mra), d_phi_x(a.basis.second.mra, 1),
    psi_t(a.basis.first), d_psi_t(a.basis.first, 1),
    psi_x(a.basis.second), d_psi_x(a.basis.second, 1),

    integral_sfsf_t(phi_t, phi_t), dd_integral_sfsf_t(d_phi_t, d_phi_t),
    integral_sfsf_x(phi_x, phi_x), dd_integral_sfsf_x(d_phi_x, d_phi_x),
    integral_ww_t(psi_t, psi_t),   dd_integral_ww_t(d_psi_t, d_psi_t),
    integral_ww_x(psi_x, psi_x),   dd_integral_ww_x(d_psi_x, d_psi_x)
{
}

template <typename T, typename Basis2D>
T
RightNormPreconditioner2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::operator()
						 (XType xtype1, int j1, int k1, XType xtype2, int j2, int k2) const
{
    T val_x, dd_val_x, val_t, dd_val_t;
    if(xtype2 == XBSpline){
        val_x = integral_sfsf_x(j2, k2, j2, k2);
        dd_val_x = dd_integral_sfsf_x(j2, k2, j2, k2);
    }
    else {
        val_x = integral_ww_x(j2, k2, j2, k2);
        dd_val_x = dd_integral_ww_x(j2, k2, j2, k2);
    }
    if(xtype1 == XBSpline){
        val_t = integral_sfsf_t(j1, k1, j1, k1);
        dd_val_t = dd_integral_sfsf_t(j1, k1, j1, k1);
    }
    else {
        val_t = integral_ww_t(j1, k1, j1, k1);
        dd_val_t = dd_integral_ww_t(j1, k1, j1, k1);
    }
    return 1./std::sqrt( (val_x+dd_val_x) + (val_t+dd_val_t)*pow2i<T>(-2*j2));
}

template <typename T, typename Basis2D>
T
RightNormPreconditioner2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::operator()(const Index2D &index) const
{
    return RightNormPreconditioner2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::operator()
			(index.index1.xtype, index.index1.j, index.index1.k, index.index2.xtype, index.index2.j, index.index2.k);
}


template <typename T, typename Basis2D, typename CGMYOperator>
RightNormPreconditioner2D<T,Basis2D,SpaceTimeCGMYOperator1D<T,Basis2D,CGMYOperator> >::RightNormPreconditioner2D
	(const SpaceTimeCGMYOperator1D<T,Basis2D,CGMYOperator> &_a)
    : a(_a),
    phi_t(a.basis.first.mra), d_phi_t(a.basis.first.mra, 1),
    phi_x(a.basis.second.mra), d_phi_x(a.basis.second.mra, 1),
    psi_t(a.basis.first), d_psi_t(a.basis.first, 1),
    psi_x(a.basis.second), d_psi_x(a.basis.second, 1),

    integral_sfsf_t(phi_t, phi_t), dd_integral_sfsf_t(d_phi_t, d_phi_t),
    integral_sfsf_x(phi_x, phi_x), dd_integral_sfsf_x(d_phi_x, d_phi_x),
    integral_ww_t(psi_t, psi_t),   dd_integral_ww_t(d_psi_t, d_psi_t),
    integral_ww_x(psi_x, psi_x),   dd_integral_ww_x(d_psi_x, d_psi_x)
{
}

template <typename T, typename Basis2D, typename CGMYOperator>
T
RightNormPreconditioner2D<T,Basis2D,SpaceTimeCGMYOperator1D<T,Basis2D,CGMYOperator> >::operator()
							 (XType xtype1, int j1, int k1, XType xtype2, int j2, int k2) const
{
    T val_x, dd_val_x, val_t, dd_val_t;
    if(xtype2 == XBSpline){
        val_x = integral_sfsf_x(j2, k2, j2, k2);
        dd_val_x = dd_integral_sfsf_x(j2, k2, j2, k2);
    }
    else {
        val_x = integral_ww_x(j2, k2, j2, k2);
        dd_val_x = dd_integral_ww_x(j2, k2, j2, k2);
    }
    if(xtype1 == XBSpline){
        val_t = integral_sfsf_t(j1, k1, j1, k1);
        dd_val_t = dd_integral_sfsf_t(j1, k1, j1, k1);
    }
	else {
		val_t = integral_ww_t(j1, k1, j1, k1);
	    dd_val_t = dd_integral_ww_t(j1, k1, j1, k1);
	}
    if (a.cgmy_x.diffusion!=0) {
    	T diffusion = a.cgmy_x.diffusion;
    	return 1./std::sqrt( (val_x+ diffusion*dd_val_x) + (val_t+dd_val_t)*pow2i<T>(-2*j2));
    }
    else {
    	T Y = a.cgmy_x.cgmy.Y;
    	return 1./std::sqrt( (val_x+std::pow(2.,Y*j2) ) + (val_t+dd_val_t)*std::pow(2.,-Y*j2));
    }

}

template <typename T, typename Basis2D, typename CGMYOperator>
T
RightNormPreconditioner2D<T,Basis2D,SpaceTimeCGMYOperator1D<T,Basis2D,CGMYOperator> >::operator()(const Index2D &index) const
{
    return RightNormPreconditioner2D<T,Basis2D,SpaceTimeCGMYOperator1D<T,Basis2D,CGMYOperator> >::operator()
			(index.index1.xtype, index.index1.j, index.index1.k, index.index2.xtype, index.index2.j, index.index2.k);
}


template <typename T, typename Basis2D>
LeftNormPreconditioner2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::LeftNormPreconditioner2D
	(const SpaceTimeHeatOperator1D<T,Basis2D> &_a)
    : a(_a),
    phi_x(a.basis.second.mra), d_phi_x(a.basis.second.mra, 1),
    psi_x(a.basis.second), d_psi_x(a.basis.second, 1),

    integral_sfsf_x(phi_x, phi_x), dd_integral_sfsf_x(d_phi_x, d_phi_x),
    integral_ww_x(psi_x, psi_x), dd_integral_ww_x(d_psi_x, d_psi_x)
{
}

template <typename T, typename Basis2D>
T
LeftNormPreconditioner2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::operator()
						(XType xtype1, int j1, int k1, XType xtype2, int j2, int k2) const
{
    T val_x, dd_val_x;
    if(xtype2 == XBSpline){
        val_x = integral_sfsf_x(j2, k2, j2, k2);
        dd_val_x = dd_integral_sfsf_x(j2, k2, j2, k2);
    }
    else {
        val_x = integral_ww_x(j2, k2, j2, k2);
        dd_val_x = dd_integral_ww_x(j2, k2, j2, k2);
    }
    return 1./std::sqrt(val_x+dd_val_x);
}

template <typename T, typename Basis2D>
T
LeftNormPreconditioner2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::operator()(const Index2D &index) const
{
    return LeftNormPreconditioner2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::operator()
			(index.index1.xtype, index.index1.j, index.index1.k, index.index2.xtype, index.index2.j, index.index2.k);
}


template <typename T, typename Basis2D, typename CGMYOperator>
LeftNormPreconditioner2D<T,Basis2D,SpaceTimeCGMYOperator1D<T,Basis2D,CGMYOperator> >::LeftNormPreconditioner2D
	(const SpaceTimeCGMYOperator1D<T,Basis2D,CGMYOperator> &_a)
    : a(_a),
    phi_x(a.basis.second.mra), d_phi_x(a.basis.second.mra, 1),
    psi_x(a.basis.second), d_psi_x(a.basis.second, 1),

    integral_sfsf_x(phi_x, phi_x), dd_integral_sfsf_x(d_phi_x, d_phi_x),
    integral_ww_x(psi_x, psi_x), dd_integral_ww_x(d_psi_x, d_psi_x)
{
}

template <typename T, typename Basis2D, typename CGMYOperator>
T
LeftNormPreconditioner2D<T,Basis2D,SpaceTimeCGMYOperator1D<T,Basis2D,CGMYOperator> >::operator()
							(XType xtype1, int j1, int k1, XType xtype2, int j2, int k2) const
{
	T val_x, dd_val_x;
	if(xtype2 == XBSpline){
		val_x = integral_sfsf_x(j2, k2, j2, k2);
	    dd_val_x = dd_integral_sfsf_x(j2, k2, j2, k2);
	}
	else {
	    val_x = integral_ww_x(j2, k2, j2, k2);
	    dd_val_x = dd_integral_ww_x(j2, k2, j2, k2);
	}
	if (a.cgmy_x.diffusion!=0) {
		T diffusion = a.cgmy_x.diffusion;
		return 1./std::sqrt(val_x+diffusion*dd_val_x);
	}
	else {
		T Y = a.cgmy_x.cgmy.Y;
		return 1./std::sqrt(val_x+std::pow(2.,Y*j2));
	}
}

template <typename T, typename Basis2D, typename CGMYOperator>
T
LeftNormPreconditioner2D<T,Basis2D,SpaceTimeCGMYOperator1D<T,Basis2D,CGMYOperator> >::operator()(const Index2D &index) const
{
    return LeftNormPreconditioner2D<T,Basis2D,SpaceTimeCGMYOperator1D<T,Basis2D,CGMYOperator> >::operator()
			(index.index1.xtype, index.index1.j, index.index1.k, index.index2.xtype, index.index2.j, index.index2.k);
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


} // namespace lawa
