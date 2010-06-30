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
#include <lawa/math/const.h>

namespace lawa {

template <QuadratureType Quad>
int
QuadratureParam<Quad>::numEvaluationsPerPiece = 100;

//--- Gauss-Legendre Quadrature-------------------------------------------------

template <typename T, typename Integrand>
Quadrature<T,Gauss,Integrand>::Quadrature(const Integrand &_integrand, int order)
    : integrand(_integrand), _order(order)
{
    _legendre();
}

template <typename T, typename Integrand>
const T
Quadrature<T,Gauss,Integrand>::operator()(T a, T b) const
{
    T ret = 0.0;
    for (int i=1; i<=_order; ++i) {
        ret += _weights(_order,i) 
               * integrand.integrand(0.5*(b-a)*_knots(_order,i)+0.5*(b+a));
    }
    ret *= 0.5*(b-a);

    return ret;
}

template <typename T, typename Integrand>
void
Quadrature<T,Gauss,Integrand>::setOrder(int order)
{
    assert(order>0);
    if (_order!=order) {
        _order = order;
        _legendre();
    }
}

template <typename T, typename Integrand>
void
Quadrature<T,Gauss,Integrand>::_legendre()
{
    T eps = Const<T>::EQUALITY_EPS;
    _knots.engine().resize(_order, _order);
    _weights.engine().resize(_order, _order);

    T x1 = -1,
      x2 =  1;
    for (int k=1; k<=_order; ++k) {
        int     m = (k+1)/2;
        T xm = 0.5 * (x2+x1),
          xl = 0.5 * (x2-x1);
        for (int i=1; i<=m; ++i) {
            T z = cos(M_PI*(i-0.25)/(k+0.5)),
              z1, pp;
            do {
                T p1 = 1.0,
                  p2 = 2.0;
                for (int j=1; j<=k; ++j) {
                    T p3 = p2;
                    p2 = p1;
                    p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
                }
                pp = k * (z*p1-p2)/(z*z-1.0);
                z1 = z;
                z = z1-p1/pp;
            } while (fabs(z-z1) > eps);
            _knots(k,i)     = xm - xl*z;
            _knots(k,k+1-i) = xm + xl*z;
            _weights(k,i)     = 2.0*xl/((1.0-z*z)*pp*pp);
            _weights(k,k+1-i) = _weights(k,i);
        }
    }
}

//--- Composite Trapezoidal rule -----------------------------------------------

template <typename T, typename Integrand>
Quadrature<T,CompositeTrapezoidal,Integrand>::
Quadrature(const Integrand &_integrand, int _n)
    : n(_n), integrand(_integrand)
{
}

template <typename T, typename Integrand>
const T
Quadrature<T,CompositeTrapezoidal,Integrand>::operator()(T a, T b) const
{
	assert(n>0);
	
    T h = (b-a) / n;
    T ret = .5 * h * integrand.integrand(a);
    a += h;
    for (int i=1; i<n; ++i, a+=h) {
        ret += h * integrand.integrand(a);
    }
    ret += .5 * h * integrand.integrand(b);

    return ret;
}

} // namespace lawa
