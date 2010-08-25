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


#ifndef LAWA_ADAPTIVE_RHSWITHPEAKS1D_H
#define LAWA_ADAPTIVE_RHSWITHPEAKS1D_H 1

#include <lawa/adaptive/index.h>
#include <lawa/adaptive/indexset.h>
#include <lawa/adaptive/coefficients.h>
#include <lawa/function.h>

namespace lawa {


template <typename T, typename Basis>
class RHSWithPeaks1D
{
    typedef typename Basis::BSplineType PrimalSpline;
    typedef typename Basis::WaveletType PrimalWavelet;

public:
    RHSWithPeaks1D(const Basis &basis, T (*f)(T), const DenseVector<Array<T> > &_singularPoints,
                   const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas, int order);

    T
    operator()(XType xtype, int j, int k) const;

    T
    operator()(const Index1D &lambda) const;

private:
    const Basis &basis;
    PrimalSpline phi;
    PrimalWavelet psi;
    Function<T> f;
    const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >    &deltas;
    Integral<T, Gauss, PrimalSpline,  Function<T> > integral_sff;
    Integral<T, Gauss, PrimalWavelet, Function<T> > integral_wf;

};


} // namespace lawa

#include <lawa/righthandsides/rhswithpeaks1d.tcc>

#endif // LAWA_ADAPTIVE_RHSWITHPEAKS1D_H
