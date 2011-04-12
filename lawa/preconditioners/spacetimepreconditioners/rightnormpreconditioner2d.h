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


#ifndef LAWA_PRECONDITIONERS_SPACETIMEPRECONDITIONERS_RIGHTNORMPRECONDITIONER2D_H
#define LAWA_PRECONDITIONERS_SPACETIMEPRECONDITIONERS_RIGHTNORMPRECONDITIONER2D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/integrals/integral.h>

namespace lawa {

template <typename T, typename Basis2D>
class RightNormPreconditioner2D
{
    typedef typename Basis2D::FirstBasisType::BSplineType PrimalSpline_t;
    typedef typename Basis2D::SecondBasisType::BSplineType PrimalSpline_x;
    typedef typename Basis2D::FirstBasisType::WaveletType PrimalWavelet_t;
    typedef typename Basis2D::SecondBasisType::WaveletType PrimalWavelet_x;

    public:
        RightNormPreconditioner2D(const Basis2D &_basis, T _s=2.);  //s=2: A: H^1 -> H^{-1}

        T
        operator()(XType xtype1, int j1, int k1, XType xtype2, int j2, int k2) const;

        T
        operator()(const Index2D &index) const;

    private:
        const Basis2D &basis;
        T             s;                //scaling for certain classes of integral operators
        PrimalSpline_t phi_t, d_phi_t;
        PrimalSpline_x phi_x, d_phi_x;
        PrimalWavelet_t psi_t, d_psi_t;
        PrimalWavelet_x psi_x, d_psi_x;

        Integral<T, Gauss, PrimalSpline_t, PrimalSpline_t>      integral_sfsf_t,
                                                             dd_integral_sfsf_t;
        Integral<T, Gauss, PrimalSpline_x, PrimalSpline_x>      integral_sfsf_x,
                                                             dd_integral_sfsf_x;
        Integral<T, Gauss, PrimalWavelet_t, PrimalWavelet_t>    integral_ww_t,
                                                             dd_integral_ww_t;
        Integral<T, Gauss, PrimalWavelet_x, PrimalWavelet_x>    integral_ww_x,
                                                             dd_integral_ww_x;
};

}   // namespace lawa

#include <lawa/preconditioners/spacetimepreconditioners/rightnormpreconditioner2d.tcc>


#endif // LAWA_PRECONDITIONERS_SPACETIMEPRECONDITIONERS_RIGHTNORMPRECONDITIONER2D_H
