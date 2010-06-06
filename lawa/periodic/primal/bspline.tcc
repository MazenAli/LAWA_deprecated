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
#include <cmath>
#include <lawa/flensforlawa.h>

#include <lawa/math/math.h>
#include <lawa/realline/primal/bspline.h>
#include <extensions/extensions.h>

namespace lawa {

using namespace flens;

template <typename T>
BSpline<T,Primal,Periodic,CDF>::BSpline(int _d)
    : d(_d), mu(d&1), deriv(0), polynomialOrder(d),
      l1(.5*(-d+mu)), l2(.5*(d+mu)),
      a(_bspline_mask<T>(d))
{
    assert(_d>0);
}

template <typename T>
BSpline<T,Primal,Periodic,CDF>::BSpline(int _d, int _deriv)
    : d(_d), mu(d&1), deriv(_deriv), polynomialOrder(d-deriv),
      l1(.5*(-d+mu)), l2(.5*(d+mu)),
      a(_bspline_mask<T>(d))
{
    assert(_d>0);
    assert(deriv>=0);
    assert(deriv<d);
}

template <typename T>
BSpline<T,Primal,Periodic,CDF>::~BSpline()
{
}

template <typename T>
T
BSpline<T,Primal,Periodic,CDF>::operator()(T x, int j, int k) const
{
    if (inner(x,support(j,k))) {
        T ret = T(0);
        x = pow2i<T>(j)*x-k - mu/2.;
        if (deriv==0) {
            x = fabs(x);
            for (int p=0; p<=ifloor(d/2.-x); ++p) {
                int sign = (p&1) ? -1 : 1;
                ret +=   sign * binomial(d, p) * pow(T(d/2.-x-p), d-1);
            }
            ret /= factorial(d - 1);
        } else {
            for (int p=0; p<=ifloor(d/2.-fabs(x)); ++p) {
                int sign = ( (p&1)==( (x>0)&&(deriv&1) ) ) ? 1 : -1;
                ret += sign * binomial(d, p)
                            * pow(T(d/2.-fabs(x)-p), d-1-deriv);
            }
            return ret/factorial(d-1-deriv);
        }
        return pow2ih<T>(j)*ret;
    }
    return T(0);
}

template <typename T>
Support<T>
BSpline<T,Primal,Periodic,CDF>::support(int j, int k) const
{
    return pow2i<T>(-j) * Support<T>(l1 + k, l2 + k);
}


template <typename T>
DenseVector<Array<T> >
BSpline<T,Primal,Periodic,CDF>::singularSupport(int j, int k) const
{
    return linspace(support(j,k).l1, support(j,k).l2, d+1);
}

template <typename T>
T
BSpline<T,Primal,Periodic,CDF>::tic(int j) const
{
    return pow2i<T>(-j);
}

template <typename T>
const DenseVector<Array<T> > &
BSpline<T,Primal,Periodic,CDF>::mask() const
{
    return a;
}

} // namespace lawa
