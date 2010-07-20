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

template <typename T, QuadratureType Quad, typename Integrand>
class Quadrature
{
};

template <QuadratureType Quad>
struct QuadratureParam
{
    static int numEvaluationsPerPiece;
};

template <typename T, typename Integrand>
class Quadrature<T,Gauss,Integrand>
{    
    public:
        Quadrature(const Integrand &_integrand, int order = 10);

        const T
        operator()(T a, T b) const;
        
        void
        setOrder(int order);

        const Integrand &integrand;
    private:
        void
        _legendre();

        int _order;
        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _knots;
        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _weights;
};

template <typename T, typename Integrand>
class Quadrature<T,CompositeTrapezoidal,Integrand>
{    
    public:
        Quadrature(const Integrand &_integrand, int _n = 1024);

        const T
        operator()(T a, T b) const;

        int n;

        const Integrand &integrand;
    private:
};

} // namespace lawa

#include <lawa/quadrature.tcc>

#endif // LAWA_QUADRATURE_H
