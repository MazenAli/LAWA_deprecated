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


#ifndef LAWA_OPERATORS_HELMHOLTZOPERATOR1D_H
#define LAWA_OPERATORS_HELMHOLTZOPERATOR1D_H 1

#include <lawa/adaptive/index.h>
#include <lawa/integrals.h>
#include <lawa/enum.h>

namespace lawa {
    
/* HELMHOLTZ OPERATOR 
 *
 *    a(u,v) =  Integral(u_x * v_x) + c * Integral(u * v)
 *
 */
template <typename T, typename Basis>
class HelmholtzOperator1D{
    
    public:
        
        const Basis& basis;
        const T c;
        
        typedef typename Basis::BSplineType PrimalSpline;
        typedef typename Basis::WaveletType PrimalWavelet;
        
        PrimalSpline phi, d_phi;
        PrimalWavelet psi, d_psi;
        
        Integral<T, Gauss, PrimalSpline, PrimalSpline> integral_sfsf, dd_integral_sfsf;
        Integral<T, Gauss, PrimalSpline, PrimalWavelet> integral_sfw, dd_integral_sfw;
        Integral<T, Gauss, PrimalWavelet, PrimalSpline> integral_wsf, dd_integral_wsf;
        Integral<T, Gauss, PrimalWavelet, PrimalWavelet> integral_ww, dd_integral_ww;
            
    public:
        HelmholtzOperator1D(const Basis& _basis, const T _c);
        HelmholtzOperator1D(const HelmholtzOperator1D<T,Basis> &a);

        T getc() const;
        const Basis& getBasis() const;
    
        T
        operator()(XType xtype1, int j1, int k1, 
                   XType xtype2, int j2, int k2) const;

        T
        operator()(const Index1D &row_index, const Index1D &col_index) const;
    
        T
        operator()(T time, 
                   XType xtype1, int j1, int k1, 
                   XType xtype2, int j2, int k2) const
        {
            return operator()(xtype1, j1, k1, xtype2, j2, k2);
        }
        
        T
        operator()(T time, const Index1D &row_index, const Index1D &col_index) const
        {
            return operator()(row_index, col_index);
        }
};
    
    
} // namespace lawa

#include <lawa/operators/helmholtzoperator1d.tcc>

#endif // LAWA_OPERATORS_HELMHOLTZOPERATOR1D_H
