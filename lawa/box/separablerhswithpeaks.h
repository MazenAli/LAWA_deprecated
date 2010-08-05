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

#ifndef LAWA_BOX_SEPARABLERHSWITHPEAKS_H
#define LAWA_BOX_SEPARABLERHSWITHPEAKS_H 1

#include <lawa/box/separablefunction.h>
#include <lawa/integrals.h>

namespace lawa {
    
template<typename T, typename Basis>
class SeparableRHSwithPeaks
{
    private:
        const Basis& basis;
        const SeparableFunction<T>& F; 
        const flens::DenseVector<flens::Array<T> > peaks_x;
        const flens::DenseVector<flens::Array<T> > peaksigns_x;
        const flens::DenseVector<flens::Array<T> > peaks_y;
        const flens::DenseVector<flens::Array<T> > peaksigns_y;
        const SeparableFunction<T>& peakF; 
        
        typedef typename Basis::FirstBasisType::BSplineType PrimalSpline_x;
        typedef typename Basis::SecondBasisType::BSplineType PrimalSpline_y;
        typedef typename Basis::FirstBasisType::WaveletType PrimalWavelet_x;
        typedef typename Basis::SecondBasisType::WaveletType PrimalWavelet_y;
        
        PrimalSpline_x phi_x;
        PrimalSpline_y phi_y;
        PrimalWavelet_x psi_x;
        PrimalWavelet_y psi_y;
        
        Integral<T, CompositeTrapezoidal, PrimalSpline_x, Function<T> > integral_sff_x, peak_integral_sff_x;
        Integral<T, CompositeTrapezoidal, PrimalSpline_y, Function<T> > integral_sff_y, peak_integral_sff_y;
        Integral<T, CompositeTrapezoidal, PrimalWavelet_x, Function<T> > integral_wf_x, peak_integral_wf_x;
        Integral<T, CompositeTrapezoidal, PrimalWavelet_y, Function<T> > integral_wf_y, peak_integral_wf_y;
                
    public:
        SeparableRHSwithPeaks(const Basis& _basis, const SeparableFunction<T>& _F, 
            const flens::DenseVector<flens::Array<T> > _peaks_x,
            const flens::DenseVector<flens::Array<T> > _peaksigns_x,
            const flens::DenseVector<flens::Array<T> > _peaks_y,
            const flens::DenseVector<flens::Array<T> > _peaksigns_y,
            const SeparableFunction<T>& _peakF);
        
        T
        operator()(bool XisSpline, int j_x, int k_x,
                   bool YisSpline, int j_y, int k_y) const;
                           
};   
    
} // namespace lawa

#include <lawa/box/separablerhswithpeaks.tcc>

#endif // LAWA_BOX_SEPARABLERHSWITHPEAKS_H