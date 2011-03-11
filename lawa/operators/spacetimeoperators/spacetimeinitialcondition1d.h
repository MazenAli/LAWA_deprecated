/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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


#ifndef LAWA_OPERATORS_SPACETIMEINITIALCONDITION1D_H
#define LAWA_OPERATORS_SPACETIMEINITIALCONDITION1D_H 1


#include <lawa/integrals.h>

namespace lawa {

template <typename T, typename Basis>
class SpaceTimeInitialCondition1D{

    public:

        const Basis& basis;

        typedef typename Basis::FirstBasisType Basis_t;
        typedef typename Basis::SecondBasisType Basis_x;

        typedef IdentityOperator1D<T, typename Basis::SecondBasisType> Reaction_x;

        Reaction_x  id_x;

    private:
        typedef typename Basis::FirstBasisType::BSplineType PrimalSpline_t;
        typedef typename Basis::SecondBasisType::BSplineType PrimalSpline_x;
        typedef typename Basis::FirstBasisType::WaveletType PrimalWavelet_t;
        typedef typename Basis::SecondBasisType::WaveletType PrimalWavelet_x;

        PrimalSpline_t phi_t;
        PrimalSpline_x phi_x;
        PrimalWavelet_t psi_t;
        PrimalWavelet_x psi_x;

        Integral<T, Gauss, PrimalSpline_x, PrimalSpline_x>   integral_sfsf_x;
        Integral<T, Gauss, PrimalSpline_x, PrimalWavelet_x>  integral_sfw_x;
        Integral<T, Gauss, PrimalWavelet_x, PrimalSpline_x>  integral_wsf_x;
        Integral<T, Gauss, PrimalWavelet_x, PrimalWavelet_x> integral_ww_x;

    public:
        SpaceTimeInitialCondition1D(const Basis& _basis);

        const Basis& getBasis() const;

        T                                                           // returns a(v,u)
        operator()(XType row_xtype_x, int j1_x, int k1_x,
                   XType col_xtype_t, int j2_t, int k2_t,
                   XType col_xtype_x, int j2_x, int k2_x) const;

        T
        operator()(const Index1D &row_index, const Index2D &col_index) const;
};


} // namespace lawa

#include <lawa/operators/spacetimeoperators/spacetimeinitialcondition1d.tcc>

#endif // LAWA_OPERATORS_SPACETIMEINITIALCONDITION1D_H
