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


#ifndef LAWA_OPERATORS_WEAKLAPLACEOPERATOR1D_H
#define LAWA_OPERATORS_WEAKLAPLACEOPERATOR1D_H 1

#include <lawa/enum.h>
#include <lawa/adaptive/index.h>
#include <lawa/integrals.h>

namespace lawa {

/* Weak Laplace OPERATOR
 *
 *    a(u,v) =  Integral(u_x * v_x)
 *
 */
template <typename T, typename Basis>
class WeakLaplaceOperator1D{

    private:

        const Basis& basis;

        typedef typename Basis::BSplineType PrimalSpline;
        typedef typename Basis::WaveletType PrimalWavelet;

        PrimalSpline d_phi;
        PrimalWavelet d_psi;

        Integral<T, Gauss, PrimalSpline, PrimalSpline>   dd_integral_sfsf;
        Integral<T, Gauss, PrimalSpline, PrimalWavelet>  dd_integral_sfw;
        Integral<T, Gauss, PrimalWavelet, PrimalSpline>  dd_integral_wsf;
        Integral<T, Gauss, PrimalWavelet, PrimalWavelet> dd_integral_ww;

    public:
        WeakLaplaceOperator1D(const Basis& _basis);

        T
        operator()(XType xtype1, int j1, int k1,
                   XType xtype2, int j2, int k2) const;

        T
        operator()(const Index1D &row_index, const Index1D &col_index) const;
};

} // namespace lawa

#include <lawa/operators/weaklaplaceoperator1d.tcc>

#endif  //LAWA_OPERATORS_WEAKLAPLACEOPERATOR1D_H 1
