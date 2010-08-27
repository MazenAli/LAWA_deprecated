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


#ifndef LAWA_OPERATORS_HELMHOLTZOPERATOR2D_H
#define LAWA_OPERATORS_HELMHOLTZOPERATOR2D_H 1

#include <lawa/enum.h>
#include <lawa/adaptive/index.h>
#include <lawa/integrals.h>

namespace lawa {

template <typename T, typename Basis2D>
class HelmholtzOperator2D{

    public:

        const Basis2D &basis;
        const T c;

        typedef typename Basis2D::FirstBasisType::BSplineType PrimalSpline_x;
        typedef typename Basis2D::SecondBasisType::BSplineType PrimalSpline_y;
        typedef typename Basis2D::FirstBasisType::WaveletType PrimalWavelet_x;
        typedef typename Basis2D::SecondBasisType::WaveletType PrimalWavelet_y;

        PrimalSpline_x phi_x, d_phi_x;
        PrimalSpline_y phi_y, d_phi_y;
        PrimalWavelet_x psi_x, d_psi_x;
        PrimalWavelet_y psi_y, d_psi_y;

        Integral<T, Gauss, PrimalSpline_x, PrimalSpline_x> integral_sfsf_x,
                                                        dd_integral_sfsf_x;
        Integral<T, Gauss, PrimalSpline_y, PrimalSpline_y> integral_sfsf_y,
                                                        dd_integral_sfsf_y;
        Integral<T, Gauss, PrimalSpline_x, PrimalWavelet_x>    integral_sfw_x,
                                                            dd_integral_sfw_x;
        Integral<T, Gauss, PrimalSpline_y, PrimalWavelet_y>    integral_sfw_y,
                                                            dd_integral_sfw_y;
        Integral<T, Gauss, PrimalWavelet_x, PrimalSpline_x>    integral_wsf_x,
                                                            dd_integral_wsf_x;
        Integral<T, Gauss, PrimalWavelet_y, PrimalSpline_y>    integral_wsf_y,
                                                            dd_integral_wsf_y;
        Integral<T, Gauss, PrimalWavelet_x, PrimalWavelet_x>    integral_ww_x,
                                                             dd_integral_ww_x;
        Integral<T, Gauss, PrimalWavelet_y, PrimalWavelet_y>    integral_ww_y,
                                                             dd_integral_ww_y;

    public:
        HelmholtzOperator2D(const Basis2D& _basis, const T _c);

        T getc() const;
        const Basis2D& getBasis() const;

        T
        operator()(XType row_xtype_x, int j1_x, int k1_x,
                   XType row_xtype_y, int j1_y, int k1_y,
                   XType col_xtype_x, int j2_x, int k2_x,
                   XType col_xtpye_y, int j2_y, int k2_y) const;

        T
        operator()(const Index2D &row_index, const Index2D &col_index) const;
};

}    // namespace lawa

#include <lawa/operators/helmholtzoperator2d.tcc>

#endif // LAWA_ADAPTIVE_HELMHOLTZOPERATOR2D_H
