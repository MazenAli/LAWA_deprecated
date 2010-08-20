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

#ifndef LAWA_RIGHTHANDSIDES_SEPARABLERHS2D_H
#define LAWA_RIGHTHANDSIDES_SEPARABLERHS2D_H 1

#include <lawa/adaptive/indexset.h>
#include <lawa/adaptive/coefficients.h>
#include <lawa/righthandsides/separablefunction2d.h>
#include <lawa/integrals.h>

namespace lawa {
    
template<typename T, typename Basis2D>
class SeparableRHS2D
{
    public:
        const Basis2D& basis;
        const SeparableFunction2D<T>& F;
        
        typedef typename Basis2D::FirstBasisType::BSplineType PrimalSpline_x;
        typedef typename Basis2D::SecondBasisType::BSplineType PrimalSpline_y;
        typedef typename Basis2D::FirstBasisType::WaveletType PrimalWavelet_x;
        typedef typename Basis2D::SecondBasisType::WaveletType PrimalWavelet_y;
        
        PrimalSpline_x phi_x;
        PrimalSpline_y phi_y;
        PrimalWavelet_x psi_x;
        PrimalWavelet_y psi_y;
        
        Integral<T, Gauss, PrimalSpline_x, Function<T> >  integral_sff_x;
        Integral<T, Gauss, PrimalSpline_y, Function<T> >  integral_sff_y;
        Integral<T, Gauss, PrimalWavelet_x, Function<T> > integral_wf_x;
        Integral<T, Gauss, PrimalWavelet_y, Function<T> > integral_wf_y;
                
    public:
        SeparableRHS2D(const Basis2D& _basis, const SeparableFunction2D<T>& _F, int order);
        
        T
        operator()(XType xtype_x, int j_x, int k_x,
                   XType xtpye_y, int j_y, int k_y) const;

        T
        operator()(const Index2D &index) const;
                           
};

template <typename T, typename RHS2D>
class SumOfRHS2D
{
private:
	const RHS2D &rhs1;
	const RHS2D &rhs2;

public:
	SumOfRHS2D(const RHS2D &rhs1, const RHS2D &rhs2);

    T
    operator()(XType xtype_x, int j_x, int k_x,
               XType xtpye_y, int j_y, int k_y) const;

    T
    operator()(const Index2D &index) const;

    Coefficients<Lexicographical,T,Index2D>
    operator()(const IndexSet<Index2D> &Lambda) const;
};
    
} // namespace lawa

#include <lawa/righthandsides/separablerhs2d.tcc>

#endif // LAWA_RIGHTHANDSIDES_SEPARABLERHS2D_H
