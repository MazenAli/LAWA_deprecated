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

#ifndef LAWA_INTERVAL_PRIMAL_BSPLINE_H
#define LAWA_INTERVAL_PRIMAL_BSPLINE_H 1

#include <flens/flens.h>

#include <lawa/bspline.h>
#include <lawa/enum.h>
#include <lawa/support.h>
#include <lawa/mra.h>

namespace lawa {

template <typename T, Construction Cons>
struct BSpline<T,Primal,Interval,Cons>
{
    typedef T ElementType;

    BSpline(const MRA<T,Primal,Interval,Cons> &_mra);

    BSpline(const MRA<T,Primal,Interval,Cons> &_mra, int _deriv);

    T
    operator()(T x, int j, int k) const;

    Support<T>
    support(int j, int k) const;

    DenseVector<Array<T> >
    singularSupport(int j, int k) const;
    
    T
    tic(int j) const;

    const MRA<T,Primal,Interval,Cons> &mra;
    const int deriv, polynomialOrder;
};

} // namespace lawa

#include <lawa/interval/primal/bspline.tcc>

#endif // LAWA_INTERVAL_PRIMAL_BSPLINE_H
