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


#ifndef LAWA_PRECONDITIONERS_PRECONDITIONERS1D_H1NORMPRECONDITIONER1D_H
#define LAWA_PRECONDITIONERS_PRECONDITIONERS1D_H1NORMPRECONDITIONER1D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/integrals/integral.h>

namespace lawa {

template <typename T, typename Basis>
class H1NormPreconditioner1D
{
    typedef typename Basis::BSplineType PrimalSpline;
    typedef typename Basis::WaveletType PrimalWavelet;

    public:
        H1NormPreconditioner1D(const Basis &_basis);

        T
        operator()(XType xtype1, int j1, int k1) const;

        T
        operator()(const Index1D &index) const;

    private:
        const Basis &basis;

        PrimalSpline phi, d_phi;
        PrimalWavelet psi, d_psi;

        Integral<T, Gauss, PrimalSpline, PrimalSpline> integral_sfsf, dd_integral_sfsf;
        Integral<T, Gauss, PrimalWavelet, PrimalWavelet> integral_ww, dd_integral_ww;

};

}   // namespace lawa

#include <lawa/preconditioners/preconditioners1d/H1normpreconditioner1d.tcc>


#endif // LAWA_PRECONDITIONERS_PRECONDITIONERS1D_H1NORMPRECONDITIONER1D_H
