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
#include <lawa/math/math.h>

namespace lawa {

//--- Gauss-Legendre Quadrature-------------------------------------------------

template <typename Integral>
Quadrature<Gauss,Integral>::Quadrature(const Integral &_integral)
    : integral(_integral), _order(-1)
{
    _legendre(_precalculatedOrder);
    setOrder(4);
}

template <typename Integral>
const typename Integral::T
Quadrature<Gauss,Integral>::operator()(T a, T b) const
{
    T ret = 0.0;
    for (int i=1; i<=_order; ++i) {
        ret += _weights(_order,i) * integral.integrand(0.5*(b-a)*_knots(_order,i)+0.5*(b+a));
    }
    ret *= 0.5*(b-a);

    return ret;
}

template <typename Integral>
void
Quadrature<Gauss,Integral>::setOrder(int order)
{
    assert(order>0);

    if (order>=_precalculatedOrder) {
        _legendre(order);
        _precalculatedOrder = order;
    }
    _order = order;
}

template <typename Integral>
int
Quadrature<Gauss,Integral>::order() const
{
    return _order;
}

template <typename Integral>
void
Quadrature<Gauss,Integral>::_legendre(int order)
{
    T eps = Const<T>::EQUALITY_EPS;
    _knots.engine().resize(order, order);
    _weights.engine().resize(order, order);

    T x1 = -1,
      x2 =  1;
    for (int k=1; k<=order; ++k) {
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

template <typename Integral>
flens::GeMatrix<flens::FullStorage<typename Integral::T,cxxblas::ColMajor> >
Quadrature<Gauss,Integral>::_weights;

template <typename Integral>
flens::GeMatrix<flens::FullStorage<typename Integral::T,cxxblas::ColMajor> >
Quadrature<Gauss,Integral>::_knots;

template <typename Integral>
int
Quadrature<Gauss,Integral>::_precalculatedOrder = 10;

//---  Trapezoidal rule -------------------------------------------------------
template <typename Integral>
Quadrature<Trapezoidal,Integral>::Quadrature(const Integral &_integral)
    : integral(_integral), _n(-1)
{
    setN(100);
}

template <typename Integral>
const typename Integral::T
Quadrature<Trapezoidal,Integral>::operator()(T a, T b) const
{
    T h = (b-a) / _n;
    T ret = .5 * h * integral.integrand(a);
    a += h;
    for (int i=1; i<_n; ++i, a+=h) {
        ret += h * integral.integrand(a);
    }
    ret += .5 * h * integral.integrand(b);

    return ret;
}

template <typename Integral>
int
Quadrature<Trapezoidal,Integral>::n() const
{
    return _n;
}

template <typename Integral>
void
Quadrature<Trapezoidal,Integral>::setN(const int n)
{
    assert(n>0);

    _n = n;
}

//---  ExpWeighted Rule -------------------------------------------------------
template <typename Integral>
Quadrature<ExpWeighted,Integral>::Quadrature(const Integral &_integral)
    : integral(_integral), _eta(-0.5*std::log(integral.function(1.))),
      _max_polynomialorder(2)
{
    if (integral._f2) {
        _max_polynomialorder = integral.first.d-1 + integral.second.d-1;
    }
    else {
        _max_polynomialorder = integral.first.d;
    }
    std::cout << "_eta = " << _eta << std::endl;
}

template <typename Integral>
const typename Integral::T
Quadrature<ExpWeighted,Integral>::operator()(T a, T b) const
{
    if (a>=b) return 0.;
    if (_max_polynomialorder<=2) {
        T x1 = a+0.25*(b-a);
        T x2 = a+0.75*(b-a);

        T f1_x1 = integral.first.generator(integral.e1)
                 ((x1-integral.left)/(integral.RightmLeft),integral.j1,integral.k1,integral.deriv1)
                 /(integral.SqrtRightmLeft);
        T f1_x2 = integral.first.generator(integral.e1)
                 ((x2-integral.left)/(integral.RightmLeft),integral.j1,integral.k1,integral.deriv1)
                 /(integral.SqrtRightmLeft);
        T f2_x1 = integral.second.generator(integral.e2)
                 ((x1-integral.left)/(integral.RightmLeft),integral.j2,integral.k2,integral.deriv2)
                 /(integral.SqrtRightmLeft);
        T f2_x2 = integral.second.generator(integral.e2)
                 ((x2-integral.left)/(integral.RightmLeft),integral.j2,integral.k2,integral.deriv2)
                 /(integral.SqrtRightmLeft);

        T alpha1 = 0.;//(f1_x2-f1_x1)/(x2-x1);
        T beta1 = f1_x1;//(x2*f1_x1-x1*f1_x2)/(x2-x1);
        T alpha2 = 0.;(f2_x2-f2_x1)/(x2-x1);
        T beta2 = f2_x1;//(x2*f2_x1-x1*f2_x2)/(x2-x1);

        //std::cout << f1_x1 << std::endl;
        //std::cout << "alpha1=" << alpha1 << ", beta1=" << beta1 << std::endl;
        //std::cout << "alpha2=" << alpha2 << ", beta2=" << beta2 << std::endl;


        T c2 = alpha1*alpha2;
        T c1 = alpha1*beta2+alpha2*beta1;
        T c0 = beta1*beta2;

        T c2_div_eta = c2/_eta;
        T One_div_2eta = 0.5*_eta;

        T ret = 0.;
        if (b<=0) {
            T Exp_2eta_a = std::exp(2*_eta*a);
            T Exp_2eta_b = std::exp(2*_eta*b);
            ret = 1./(2*_eta)*c0*(Exp_2eta_b-Exp_2eta_a);
            /*
            ret +=   c2/(2*_eta)*(b*b*Exp_2eta_b-a*a*Exp_2eta_a)
                   + 1./(2*_eta)*(c1-c2_div_eta)*(b*Exp_2eta_b-a*Exp_2eta_a)
                   + 1./(2*_eta)*(c0-(c1-c2_div_eta)/(2*_eta))*(Exp_2eta_b-Exp_2eta_a);
             */
        }
        else if (a<0 && b>0) {
            return -100.;
            T Exp_2eta_a = std::exp(2*_eta*a);
            T Exp_m_2eta_b = std::exp(-2*_eta*b);
            ret += - c2/(2*_eta)*(b*b*Exp_m_2eta_b+a*a*Exp_2eta_a)
                   - 1./(2*_eta)*(c1-c2_div_eta)*a*Exp_2eta_a
                   - 1./(2*_eta)*(c1+c2_div_eta)*b*Exp_m_2eta_b
                   + 1./(2*_eta)*(c0-(c1-c2_div_eta)/(2*_eta))*(1-Exp_2eta_a)
                   + 1./(2*_eta)*(c0+(c1+c2_div_eta)/(2*_eta))*(1-Exp_m_2eta_b);

        }
        else if (a>=0) {
            T Exp_m_2eta_a = std::exp(-2*_eta*a);
            T Exp_m_2eta_b = std::exp(-2*_eta*b);
            ret = 1./(2*_eta)*c0*(Exp_m_2eta_a-Exp_m_2eta_b);
            /*
            ret +=   c2/(2*_eta)*(a*a*Exp_m_2eta_a-b*b*Exp_m_2eta_b)
                   + 1./(2*_eta)*(c1+c2_div_eta)*(a*Exp_m_2eta_a-b*Exp_m_2eta_b)
                   + 1./(2*_eta)*(c0+(c1+c2_div_eta)/(2*_eta))*(Exp_m_2eta_a-Exp_m_2eta_b);
            */
        }
        return ret;
    }
    else {
        assert(0);
        return 0.;
    }
}
} // namespace lawa

