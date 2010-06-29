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

namespace lawa {

using namespace flens;

template <typename T>
Wavelet<T,Primal,R,CDF>::Wavelet(int _d, int _d_)
    : d(_d), d_(_d_), mu(d&1), l1((2-d-d_)/2), l2((d+d_)/2),    
      deriv(0), polynomialOrder(_d),
      vanishingMoments(_d_), b(mask(d,d_)),
      phi(d), phi_(d,d_)
{
    assert(d<=d_);
    assert(((d+d_)&1)==0);
}

template <typename T>
Wavelet<T,Primal,R,CDF>::Wavelet(int _d, int _d_, int _deriv)
    : d(_d), d_(_d_), mu(d&1), l1((2-d-d_)/2), l2((d+d_)/2),    
      deriv(_deriv), polynomialOrder(_d-deriv),
      vanishingMoments(_d_), b(mask(d,d_)),
      phi(d), phi_(d,d_)
{
    assert(d<=d_);
    assert(((d+d_)&1)==0);
    assert(deriv>=0);
}

template <typename T>
Wavelet<T,Primal,R,CDF>::Wavelet(const BSpline<T,Primal,R,CDF> &_phi,    
                                 const BSpline<T,Dual,R,CDF> &_phi_)
    : d(_phi.d), d_(_phi_.d_), mu(d&1), l1((2-d-d_)/2), l2((d+d_)/2),    
      deriv(0), polynomialOrder(d),
      vanishingMoments(d_), b(mask(d,d_)),
      phi(_phi), phi_(_phi_)
      
{
}

template <typename T>
Wavelet<T,Primal,R,CDF>::Wavelet(const BSpline<T,Primal,R,CDF> &_phi,    
                                 const BSpline<T,Dual,R,CDF> &_phi_,
                                 int _deriv)
    : d(_phi.d), d_(_phi_.d_), mu(d&1), l1((2-d-d_)/2), l2((d+d_)/2),    
      deriv(_deriv), polynomialOrder(d-deriv),
      vanishingMoments(d_), b(mask(d,d_)),
      phi(_phi), phi_(_phi_)
{    
    assert(deriv>=0);
}

template <typename T>
T
Wavelet<T,Primal,R,CDF>::operator()(T x, int j, int k) const
{
    T ret = T(0);
    x = pow2i<T>(j)*x-k;
    for (int i=b.firstIndex(); i<=b.lastIndex(); ++i) {
        ret += b(i)*phi(2*x-i, 0, 0);
    }
    return pow2i<T>(deriv*(j+1)) * pow2ih<T>(j) * ret;
}

template <typename T>
Support<T>
Wavelet<T,Primal,R,CDF>::support(int j, int k) const
{
    return pow2i<T>(-j) * Support<T>(l1+k, l2+k);
}

template <typename T>
DenseVector<Array<T> >
Wavelet<T,Primal,R,CDF>::singularSupport(int j, int k) const
{
    return linspace(support(j,k).l1, support(j,k).l2, 2*(d+d_)-1);
}

template <typename T>
T
Wavelet<T,Primal,R,CDF>::tic(int j) const
{
    return pow2i<T>(-(j+1));
}

template <typename T>
const DenseVector<Array<T> > &
Wavelet<T,Primal,R,CDF>::mask() const
{
    return b;
}

template <typename T>
DenseVector<Array<T> >
Wavelet<T,Primal,R,CDF>::mask(int d, int d_)
{
    assert(d<=d_);
    assert(((d+d_)&1)==0);

    int mu = d & 1;
    BSpline<T,Dual,R,CDF> phi_(d,d_);
    DenseVector<Array<T> > b(_(2-(d+mu)/2-d_, (d-mu)/2+d_));
    for (int k=b.firstIndex(); k<=b.lastIndex(); ++k) {
        int sign = (k&1) ? -1 : 1;
        b(k) = sign * phi_.a_(1-k);
    }
    return b;
}

} // namespace lawa
