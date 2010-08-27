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

namespace lawa {

template<typename T, typename Basis2D>
SeparableRHS2DwithPeaks<T, Basis2D>::SeparableRHS2DwithPeaks(const Basis2D& _basis, const SeparableFunction2D<T>& _F,
        const flens::DenseVector<flens::Array<T> > _peaks_x,
        const flens::DenseVector<flens::Array<T> > _peaksigns_x,
        const flens::DenseVector<flens::Array<T> > _peaks_y,
        const flens::DenseVector<flens::Array<T> > _peaksigns_y,
        const SeparableFunction2D<T>& _peakF)
    : basis(_basis), F(_F), peaks_x(_peaks_x), peaksigns_x(_peaksigns_x),
      peaks_y(_peaks_y), peaksigns_y(_peaksigns_y), peakF(_peakF), phi_x(_basis.first.mra), 
      phi_y(_basis.second.mra), psi_x(_basis.first), psi_y(_basis.second),
      integral_sff_x(phi_x, _F.F_x), peak_integral_sff_x(phi_x, _peakF.F_x),
      integral_sff_y(phi_y, _F.F_y), peak_integral_sff_y(phi_y, _peakF.F_y),
      integral_wf_x(psi_x, _F.F_x), peak_integral_wf_x(psi_x, _peakF.F_x),
      integral_wf_y(psi_y, _F.F_y), peak_integral_wf_y(psi_y, _peakF.F_y)
{  
}

template<typename T, typename Basis2D>
T
SeparableRHS2DwithPeaks<T, Basis2D>::operator()(XType xtype_x, int j_x, int k_x,
                                                XType xtype_y, int j_y, int k_y) const
{
    T val_x = 0;
    T val_y = 0;
    T val_peaks = 0;
    if(xtype_x == XBSpline){
        val_x = integral_sff_x(j_x, k_x);
        if(xtype_y == XBSpline){
            val_y = integral_sff_y(j_y, k_y);
            for(int i = peaks_x.firstIndex(); i <= peaks_x.lastIndex(); ++i){
                val_peaks += peaksigns_x(i)* phi_x(peaks_x(i), j_x, k_x) * peak_integral_sff_y(j_y, k_y);    
            }
            for(int i = peaks_y.firstIndex(); i <= peaks_y.lastIndex(); ++i){
                val_peaks += peak_integral_sff_x(j_x, k_x) * peaksigns_y(i)* phi_y(peaks_y(i), j_y, k_y);    
            }
        }
        else{
            val_y = integral_wf_y(j_y, k_y);
            for(int i = peaks_x.firstIndex(); i <= peaks_x.lastIndex(); ++i){
                val_peaks += peaksigns_x(i)* phi_x(peaks_x(i), j_x, k_x) * peak_integral_wf_y(j_y, k_y);    
            }
            for(int i = peaks_y.firstIndex(); i <= peaks_y.lastIndex(); ++i){
                val_peaks += peak_integral_sff_x(j_x, k_x) * peaksigns_y(i)* psi_y(peaks_y(i), j_y, k_y);    
            }
        }
        
    }
    else{
        val_x = integral_wf_x(j_x, k_x);
        if(xtype_y == XBSpline){
            val_y = integral_sff_y(j_y, k_y);
            for(int i = peaks_x.firstIndex(); i <= peaks_x.lastIndex(); ++i){
                val_peaks += peaksigns_x(i)* psi_x(peaks_x(i), j_x, k_x) * peak_integral_sff_y(j_y, k_y);    
            }
            for(int i = peaks_y.firstIndex(); i <= peaks_y.lastIndex(); ++i){
                val_peaks += peak_integral_wf_x(j_x, k_x) * peaksigns_y(i)* phi_y(peaks_y(i), j_y, k_y);    
            }
        }
        else{
            val_y = integral_wf_y(j_y, k_y);
            for(int i = peaks_x.firstIndex(); i <= peaks_x.lastIndex(); ++i){
                val_peaks += peaksigns_x(i)* psi_x(peaks_x(i), j_x, k_x) * peak_integral_wf_y(j_y, k_y);    
            }
            for(int i = peaks_y.firstIndex(); i <= peaks_y.lastIndex(); ++i){
                val_peaks += peak_integral_wf_x(j_x, k_x) * peaksigns_y(i)* psi_y(peaks_y(i), j_y, k_y);    
            }
        }
    }
    
    return val_x * val_y + val_peaks;
}

template<typename T, typename Basis2D>
T
SeparableRHS2DwithPeaks<T, Basis2D>::operator()(const Index2D &index) const
{
    return SeparableRHS2DwithPeaks<T, Basis2D>::operator()(index.index1.xtype, index.index1.j, index.index1.k,
                                                           index.index2.xtype, index.index2.j, index.index2.k);
}

} // namespace lawa
