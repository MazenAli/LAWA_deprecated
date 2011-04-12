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

#ifndef LAWA_RIGHTHANDSIDES_SMOOTHRHSWITHALIGNEDSING2D_H
#define LAWA_RIGHTHANDSIDES_SMOOTHRHSWITHALIGNEDSING2D_H 1

#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/functiontypes/functionnd.h>
#include <lawa/integrals/integrals2d.h>

namespace lawa {

template<typename T, typename Basis2D, QuadratureType Quad>
class SmoothRHSWithAlignedSing2D
{
	const Basis2D& basis;

    typedef typename Basis2D::FirstBasisType::BSplineType PrimalSpline_x;
    typedef typename Basis2D::SecondBasisType::BSplineType PrimalSpline_y;
    typedef typename Basis2D::FirstBasisType::WaveletType PrimalWavelet_x;
    typedef typename Basis2D::SecondBasisType::WaveletType PrimalWavelet_y;

    PrimalSpline_x phi_x;
    PrimalSpline_y phi_y;
    PrimalWavelet_x psi_x;
    PrimalWavelet_y psi_y;

    Integral2D<T,Quad,PrimalSpline_x, PrimalSpline_y>  integral2d_phi_x_phi_y;
    Integral2D<T,Quad,PrimalWavelet_x,PrimalSpline_y>  integral2d_psi_x_phi_y;
    Integral2D<T,Quad,PrimalSpline_x, PrimalWavelet_y> integral2d_phi_x_psi_y;
    Integral2D<T,Quad,PrimalWavelet_x,PrimalWavelet_y> integral2d_psi_x_psi_y;

public:
    SmoothRHSWithAlignedSing2D(const Basis2D& _basis, const Function2D<T>& _F, int order);

    SmoothRHSWithAlignedSing2D(const Basis2D& _basis, const Function2D<T>& _F, int deriv_x,
							   int deriv_y, int order);

    T
    operator()(XType xtype_x, int j_x, int k_x,
               XType xtype_y, int j_y, int k_y) const;

    T
    operator()(const Index2D &index) const;

};


}	//namespace lawa

#include <lawa/righthandsides/smoothrhswithalignedsing2d.tcc>


#endif //LAWA_RIGHTHANDSIDES_SMOOTHRHSWITHALIGNEDSING2D_H 1
