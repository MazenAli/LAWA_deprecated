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

#ifndef LAWA_RIGHTHANDSIDES_SEPARABLERHS2D_H
#define LAWA_RIGHTHANDSIDES_SEPARABLERHS2D_H 1

#include <lawa/functiontypes/separablefunction2d.h>
#include <lawa/integrals/integrals.h>

namespace lawa {

template<typename T, typename Basis2D>
class SeparableRHSWithPeaks2D
{
    public:
        SeparableRHSWithPeaks2D(const Basis2D& _basis, const SeparableFunction2D<T>& _F,
                                const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_x,
                                const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_y,
                                int order);

        T
        operator()(XType xtype_x, int j_x, int k_x,
                   XType xtype_y, int j_y, int k_y) const;

        T
        operator()(const Index2D &index) const;


    private:
        const Basis2D& basis;
        const SeparableFunction2D<T>& F;
        const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_x;
        const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_y;

        IntegralF<Gauss, Basis2D::FirstBasisType>  integralf_x;
        IntegralF<Gauss, Basis2D::SecondBasisType> integralf_y;

};

} // namespace lawa

#include <lawa/righthandsides/separablerhswithpeaks2d.tcc>

#endif // LAWA_RIGHTHANDSIDES_SEPARABLERHS2D_H
