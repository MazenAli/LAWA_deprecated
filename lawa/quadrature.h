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

#ifndef LAWA_QUADRATURE_H
#define LAWA_QUADRATURE_H 1

#include <lawa/enum.h>
#include <lawa/flensforlawa.h>

namespace lawa {

template <typename T, QuadratureType Quad, typename First, typename Second>
    class Integral;

//-----------------------------------------------------------------------------

template <typename T, QuadratureType Quad, typename First, typename Second>
class Quadrature
{
};

template <typename T, typename First, typename Second>
class Quadrature<T,Gauss,First,Second>
{
    public:
        Quadrature(const Integral<T,Gauss,First,Second> &_integral);

        const T
        operator()(T a, T b) const;

        int
        order() const;

        void
        setOrder(int order);

        const Integral<T,Gauss,First,Second> &integral;
    private:
        void
        _legendre();

        int _order;
        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _knots;
        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _weights;
};

template <typename First, typename Second>
    int
    defaultGaussOrder(const First &first, const Second &second);

//-----------------------------------------------------------------------------

template <typename T, typename First, typename Second>
class Quadrature<T,Trapezoidal,First,Second>
{
    public:
        Quadrature(const Integral<T,Trapezoidal,First,Second> &_integral);

        const T
        operator()(T a, T b) const;

        int
        n() const;

        void
        setN(int n);

        const Integral<T,Trapezoidal,First,Second>  &integral;
    private:
        int _n;
};

template <typename First, typename Second>
    int
    defaultTrapezoidalN(const First &first, const Second &second);

} // namespace lawa

#include <lawa/quadrature.tcc>

#endif // LAWA_QUADRATURE_H
