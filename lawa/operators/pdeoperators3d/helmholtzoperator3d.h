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


#ifndef LAWA_OPERATORS_PDEOPERATORS3D_HELMHOLTZOPERATOR3D_H
#define LAWA_OPERATORS_PDEOPERATORS3D_HELMHOLTZOPERATOR3D_H 1

#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/operators/pdeoperators1d/laplaceoperator1d.h>
#include <lawa/operators/pdeoperators1d/identityoperator1d.h>
#include <lawa/integrals/integral.h>

namespace lawa {

template <typename T, typename Basis3D>
class HelmholtzOperator3D{

    public:

        const Basis3D &basis;
        const T c;

        typedef typename Basis3D::FirstBasisType  Basis_x;
        typedef typename Basis3D::SecondBasisType Basis_y;
        typedef typename Basis3D::ThirdBasisType  Basis_z;

        typedef LaplaceOperator1D<T, typename Basis3D::FirstBasisType>   Diffusion_x;
        typedef IdentityOperator1D<T, typename Basis3D::FirstBasisType> 	     Reaction_x;
        typedef LaplaceOperator1D<T, typename Basis3D::SecondBasisType>  Diffusion_y;
        typedef IdentityOperator1D<T, typename Basis3D::SecondBasisType> 	     Reaction_y;
        typedef LaplaceOperator1D<T, typename Basis3D::ThirdBasisType>   Diffusion_z;
        typedef IdentityOperator1D<T, typename Basis3D::ThirdBasisType> 	     Reaction_z;

        typedef typename Basis3D::FirstBasisType::BSplineType  PrimalSpline_x;
        typedef typename Basis3D::SecondBasisType::BSplineType PrimalSpline_y;
        typedef typename Basis3D::ThirdBasisType::BSplineType  PrimalSpline_z;
        typedef typename Basis3D::FirstBasisType::WaveletType  PrimalWavelet_x;
        typedef typename Basis3D::SecondBasisType::WaveletType PrimalWavelet_y;
        typedef typename Basis3D::ThirdBasisType::WaveletType  PrimalWavelet_z;

        Diffusion_x dd_x;
        Reaction_x  id_x;
        Diffusion_y dd_y;
        Reaction_y  id_y;
        Diffusion_y dd_z;
        Reaction_y  id_z;

        PrimalSpline_x phi_x, d_phi_x;
        PrimalSpline_y phi_y, d_phi_y;
        PrimalSpline_z phi_z, d_phi_z;
        PrimalWavelet_x psi_x, d_psi_x;
        PrimalWavelet_y psi_y, d_psi_y;
        PrimalWavelet_z psi_z, d_psi_z;

        Integral<T, Gauss, PrimalSpline_x, PrimalSpline_x> integral_sfsf_x,
                                                        dd_integral_sfsf_x;
        Integral<T, Gauss, PrimalSpline_y, PrimalSpline_y> integral_sfsf_y,
                                                        dd_integral_sfsf_y;
        Integral<T, Gauss, PrimalSpline_z, PrimalSpline_z> integral_sfsf_z,
                                                        dd_integral_sfsf_z;

        Integral<T, Gauss, PrimalSpline_x, PrimalWavelet_x>    integral_sfw_x,
                                                            dd_integral_sfw_x;
        Integral<T, Gauss, PrimalSpline_y, PrimalWavelet_y>    integral_sfw_y,
                                                            dd_integral_sfw_y;
        Integral<T, Gauss, PrimalSpline_z, PrimalWavelet_z>    integral_sfw_z,
                                                            dd_integral_sfw_z;

        Integral<T, Gauss, PrimalWavelet_x, PrimalSpline_x>    integral_wsf_x,
                                                            dd_integral_wsf_x;
        Integral<T, Gauss, PrimalWavelet_y, PrimalSpline_y>    integral_wsf_y,
                                                            dd_integral_wsf_y;
        Integral<T, Gauss, PrimalWavelet_z, PrimalSpline_z>     integral_wsf_z,
                                                             dd_integral_wsf_z;

        Integral<T, Gauss, PrimalWavelet_x, PrimalWavelet_x>    integral_ww_x,
                                                             dd_integral_ww_x;
        Integral<T, Gauss, PrimalWavelet_y, PrimalWavelet_y>    integral_ww_y,
                                                             dd_integral_ww_y;
        Integral<T, Gauss, PrimalWavelet_z, PrimalWavelet_z>    integral_ww_z,
                                                             dd_integral_ww_z;

    public:
        HelmholtzOperator3D(const Basis3D& _basis, const T _c);

        T getc() const;
        const Basis3D& getBasis() const;

        T
        operator()(XType row_xtype_x, int j1_x, int k1_x,
                   XType row_xtype_y, int j1_y, int k1_y,
                   XType row_xtype_z, int j1_z, int k1_z,
                   XType col_xtype_x, int j2_x, int k2_x,
                   XType col_xtpye_y, int j2_y, int k2_y,
                   XType col_xtpye_z, int j2_z, int k2_z) const;

        T
        operator()(const Index3D &row_index, const Index3D &col_index) const;
};

}    // namespace lawa

#include <lawa/operators/pdeoperators3d/helmholtzoperator3d.tcc>

#endif //  LAWA_OPERATORS_PDEOPERATORS3D_HELMHOLTZOPERATOR3D_H
