/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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

#ifndef LAWA_PERIODIC_DUAL_WAVELET_H
#define LAWA_PERIODIC_DUAL_WAVELET_H 1

#include <lawa/flensforlawa.h>
#include <lawa/periodic/primal/bspline.h>
#include <lawa/periodic/dual/bspline.h>
#include <lawa/wavelet.h>

namespace lawa {

using namespace flens;

template <typename T>
class Wavelet<T,Dual,Periodic,CDF>
{
    public:
        typedef T ElementType;

        Wavelet(int _d, int _d_);

        Wavelet(const BSpline<T,Primal,Periodic,CDF> &_phi,
                const BSpline<T,Dual,Periodic,CDF> &_phi_);

        T
        operator()(T x, int j, int k) const;

        Support<T>
        support(int j, int k) const;

        const DenseVector<Array<T> > &
        mask() const;
        
        static DenseVector<Array<T> >
        mask(int d, int d_);
        
        const int d, d_, mu;
        const int l1_, l2_;
        const DenseVector<Array<T> > b_;
        const BSpline<T,Primal,Periodic,CDF> phi;
        const BSpline<T,Dual,Periodic,CDF> phi_;        
};

} // namespace lawa

#include <lawa/periodic/dual/wavelet.tcc>

#endif // LAWA_PERIODIC_DUAL_WAVELET_H