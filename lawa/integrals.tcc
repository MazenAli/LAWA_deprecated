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
#include <lawa/aux/typetraits.h>
#include <lawa/enum.h>
#include <lawa/bspline.h>
#include <lawa/math/math.h>
#include <lawa/param.h>
#include <lawa/support.h>
#include <lawa/wavelet.h>
#include <lawa/periodic/integrals.h>

namespace lawa {

//--- primal * primal
template <typename T, typename First, typename Second>
typename RestrictTo<BothPrimal<First,Second>::value, T>::Type
_integrate(const Integral<T,Gauss,First,Second> &integral)
{
    const First &first = integral.first;
    const Second &second = integral.second;
    // quadrature with minimal order to guarantee exactness of integrals.

    /*static*/ Quadrature<T,Gauss,Integral<T,Gauss,First,Second> > quadrature(
                  integral,
                  iceil((first.polynomialOrder+second.polynomialOrder-1)/2.)+1);
    //quadrature.setOrder(iceil((first.polynomialOrder+second.polynomialOrder-1)/2.)+1);

    // the (minimal) width of the polynomial pieces.
    T unit = std::min(first.tic(integral.j1), second.tic(integral.j2));

    T ret = 0.;
    Support<T> common;
    if (overlap(first.support(integral.j1,integral.k1),
               second.support(integral.j2,integral.k2),common)) {
        T a = common.l1;
        for (T b=a+unit; b<=common.l2; b+=unit) {
            ret += quadrature(a,b);
            a = b;
         }
     }
     return ret;
}

//--- primal * function
template <typename T, QuadratureType Quad, typename First, typename Second>
typename RestrictTo<IsPrimal<First>::value && !PrimalOrDual<Second>::value, T>::Type
_integrate(const Integral<T,Quad,First,Second> &integral)
{
    
    const First &first = integral.first;
    const Second &second = integral.second;
    /*static*/ Quadrature<T,Quad,Integral<T,Quad,First,Second> > quadrature(integral);

    
    // merge singular points of bspline/wavelet and function to one list.
    /* -> implizite Annahme: second.singularPoints sind schon sortiert!! */
   // T unit = first.tic(integral.j1);
    DenseVector<Array<T> > firstSingularPoints 
                               = first.singularSupport(integral.j1,integral.k1);
    int nFirst = firstSingularPoints.length(),
        nSecond = second.singularPoints.length();
    DenseVector<Array<T> > singularPoints(nFirst + nSecond);

    std::merge(firstSingularPoints.engine().data(),
               firstSingularPoints.engine().data() + nFirst,
               second.singularPoints.engine().data(),
               second.singularPoints.engine().data() + nSecond,
               singularPoints.engine().data());
               
    T ret = 0.0;
    for (int i=singularPoints.firstIndex(); i<singularPoints.lastIndex(); ++i) {
        ret += quadrature(singularPoints(i),singularPoints(i+1));
    }
    return ret;
}

//--- primal * dual
template <typename T, QuadratureType Quad, typename First, typename Second>
typename RestrictTo<IsPrimal<First>::value && IsDual<Second>::value, T>::Type
_integrate(const Integral<T,Quad,First,Second> &integral)
{
    const First &first = integral.first;
    const Second &second = integral.second;
    /*static*/ Quadrature<T,Quad,Integral<T,Quad,First,Second> > quadrature(integral);
    
    Support<T> common;

    if (overlap(first.support(integral.j1,integral.k1),
                second.support(integral.j2,integral.k2),
                common)) {
        quadrature.n =  pow2i<T>(Param<Second>::resolution) * 1./first.tic(integral.j1);
        return quadrature(common.l1, common.l2);
    } else {
        return 0;
    }
}

//--- dual * dual
template <typename T, QuadratureType Quad, typename First, typename Second>
typename RestrictTo<BothDual<First,Second>::value, T>::Type
_integrate(const Integral<T,Quad,First,Second> &integral)
{
    const First &first = integral.first;
    const Second &second = integral.second;
   /*static*/ Quadrature<T,Quad,Integral<T,Quad,First,Second> > quadrature(integral);

    Support<T> common;
    if (overlap(first.support(0,integral.k1),
                second.support(0,integral.k2),
                common) > 0) {
        quadrature.n = pow2i<T>(Param<First>::resolution) * common.length();
        if (IsWavelet<First>::value || IsWavelet<Second>::value) {
            quadrature.n *= 2;
        }
        return quadrature(common.l1, common.l2);
    } else {
        return 0;
    }
}

//--- dual * function
template <typename T, QuadratureType Quad, typename First, typename Second>
typename RestrictTo<IsDual<First>::value && !PrimalOrDual<Second>::value, T>::Type
_integrate(const Integral<T,Quad,First,Second> &integral)
{
    const First &first = integral.first;
    const Second &second = integral.second;
    /*static*/ Quadrature<T,Quad,Integral<T,Quad,First,Second> > quadrature(integral);
    quadrature.n = pow2i<T>(Param<First>::resolution)*(first.mask().length()-1);

    Support<T> supp = first.support(integral.j1,integral.k1);
    DenseVector<Array<T> > firstSupport(2);
    firstSupport = supp.l1, supp.l2;
    // merge singular points of bspline/wavelet and function to one list.
    int nFirst = firstSupport.length(),
        nSecond = second.singularPoints.length();
    DenseVector<Array<T> > singularPoints(nFirst + nSecond);

    std::merge(firstSupport.engine().data(),
               firstSupport.engine().data() + nFirst,
               second.singularPoints.engine().data(),
               second.singularPoints.engine().data() + nSecond,
               singularPoints.engine().data());

    T ret = 0.0;
    for (int i=singularPoints.firstIndex(); i<singularPoints.lastIndex(); ++i) {
        ret += quadrature(singularPoints(i),singularPoints(i+1));
    }
    return ret;
}

//------------------------------------------------------------------------------

//--- primal * primal  or  dual * dual
template <typename T, QuadratureType Quad, typename First, typename Second>
typename RestrictTo<PrimalOrDual<First>::value
                 && PrimalOrDual<Second>::value, T>::Type
_integrand(const Integral<T,Quad,First,Second> &integral, T x)
{
    const First &first = integral.first;
    const Second &second = integral.second;
    return first(x,integral.j1,integral.k1) * second(x,integral.j2,integral.k2);
}

//--- primal or dual * anything else
template <typename T, QuadratureType Quad, typename First, typename Second>
typename RestrictTo<PrimalOrDual<First>::value
                && !PrimalOrDual<Second>::value, T>::Type
_integrand(const Integral<T,Quad,First,Second> &integral, T x)
{
    const First &first = integral.first;
    const Second &second = integral.second;
    return first(x,integral.j1,integral.k1) * second(x);
}

//------------------------------------------------------------------------------

template <typename T, QuadratureType Quad, typename First, typename Second>
Integral<T,Quad,First,Second>::Integral(const First &_first,
                                        const Second &_second)
    : first(_first), second(_second)
{
}

template <typename T, QuadratureType Quad, typename First, typename Second>
T
Integral<T,Quad,First,Second>::operator()(int _j1, int _k1, int _j2, int _k2) const
{
    j1 = _j1; k1 = _k1; j2 = _j2; k2 = _k2;
    return _integrate(*this);
}

template <typename T, QuadratureType Quad, typename First, typename Second>
T
Integral<T,Quad,First,Second>::operator()(int _j1, int _k1) const
{
    j1 = _j1; k1 = _k1;
    return _integrate(*this);
}

template <typename T, QuadratureType Quad, typename First, typename Second>
T
Integral<T,Quad,First,Second>::integrand(T x) const
{
    return _integrand(*this, x);
}

} // namespace lawa
