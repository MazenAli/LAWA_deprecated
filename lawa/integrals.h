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

#ifndef LAWA_INTEGRALS_H
#define LAWA_INTEGRALS_H 1

#include <lawa/quadrature.h>

namespace lawa {

template <typename T, QuadratureType Quad, typename First, typename Second>
struct Integral
{
    Integral(const First &_first, const Second &_second);

    T
    operator()(int _j1, int _k1, int _j2, int _k2);

    T
    operator()(int _j1, int _k1);

    T
    integrand(T x) const;

    const First &first;
    const Second &second;
    int j1, k1, j2, k2;
};

} // namespace lawa

#include <lawa/integrals.tcc>

#endif // LAWA_INTEGRALS_H