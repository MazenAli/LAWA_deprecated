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

#ifndef LAWA_INTERVAL_PRIMAL_WAVELET_H
#define LAWA_INTERVAL_PRIMAL_WAVELET_H 1

#include <lawa/support.h>
#include <lawa/basis.h>
#include <lawa/wavelet.h>

namespace lawa {

template <typename T, Construction Cons>
struct Wavelet<T,Primal,Interval,Cons>
{
    typedef T ElementType;

    Wavelet(const Basis<T,Primal,Interval,Cons> &_basis);

    Wavelet(const Basis<T,Primal,Interval,Cons> &_basis, int _deriv);

    T
    operator()(T x, int j, int k) const;

    Support<T>
    support(int j, int k) const;

    DenseVector<Array<T> >
    singularSupport(int j, int k) const;

    int
    vanishingMoments(int j, int k) const;
    
    T
    tic(int j) const;

    const Basis<T,Primal,Interval,Cons> &basis;
    const BSpline<T,Primal,Interval,Cons> phi;
    const int deriv, polynomialOrder;
};

} // namespace lawa

#include <lawa/interval/primal/wavelet.tcc>

#endif // LAWA_INTERVAL_PRIMAL_WAVELET_H
