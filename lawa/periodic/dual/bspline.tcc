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

#include <cassert>
#include <lawa/flensforlawa.h>
#include <lawa/math/polynomial.h>
#include <lawa/param.h>
#include <lawa/realline/subdivision.h>

namespace lawa {

using namespace flens;

template <typename T>
    DenseVector<Array<T> >
    _bspline_mask(int d, int d_);

template <typename T>
BSpline<T,Dual,Periodic,CDF>::BSpline(int _d, int _d_)
    : d(_d), d_(_d_), mu(d&1),
      l1_(.5*(-d+mu)-d_+1), l2_(.5*(d+mu)+d_-1),
      a_(_bspline_mask<T>(d,d_))
{
    assert(d>0);
    assert(d_>=d);
    assert(((d+d_)&1)==0);
}

template <typename T>
BSpline<T,Dual,Periodic,CDF>::~BSpline()
{
}

template <typename T>
T
BSpline<T,Dual,Periodic,CDF>::operator()(T x, int j, int k) const
{
    int resolution = Param<BSpline<T,Dual,Periodic,CDF> >::resolution;
    // we precompute values for dual B-spline on first call ...
    static DenseVector<Array<T> > values;
    static int storedD = 0, storedD_ = 0;
    static int storedResolution = resolution;
    Support<T> supp = support(j,k);
    if (!inner(x,supp)) {
        return 0;
    }

    // we need to recalculate for different B-spline than stored one.
    if ((d!=storedD) || (d_!=storedD_) || (storedResolution!=resolution)) {
        storedD  = d;
        storedD_ = d_;
        storedResolution = resolution;
        BSpline<T,Dual,R,CDF> phi_(d,d_);
        subdivide(phi_, resolution, values);
    }

    // 'revert to reference B-spline i.e. j=k=0
    x *= pow2i<T>(j);
    x -= k;
    x *= pow2i<T>(resolution);

    assert(x>=pow2i<T>(resolution)*l1_);
    assert(x<=pow2i<T>(resolution)*l2_);

    // use linear interpolation between neighboring grid points.
    return pow2ih<T>(j) * values(x);//(ifloor(values(x))+iceil(values(x)))/2;
}


template <typename T>
Support<T>
BSpline<T,Dual,Periodic,CDF>::support(int j, int k) const
{
    return pow2i<T>(-j) * Support<T>(l1_+k, l2_+k);
}

template <typename T>
const DenseVector<Array<T> > &
BSpline<T,Dual,Periodic,CDF>::mask() const
{
    return a_;
}

} // namespace lawa
