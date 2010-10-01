/*
 LAWA - Library for Adaptive Wavelet Applications.
 Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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

#ifndef LAWA_RIGHTHANDSIDES_INTEGRALS2D_H
#define LAWA_RIGHTHANDSIDES_INTEGRALS2D_H 1

#include <lawa/righthandsides/functionnd.h>
#include <lawa/righthandsides/quadrature2d.h>

namespace lawa {

template <typename T, QuadratureType Quad, typename BasisFunc_x, typename BasisFunc_y>
struct Integral2D
{
	const Function2D<T>& F;
	const BasisFunc_x &phi_x;
	const BasisFunc_y &phi_y;
    mutable int j_x, k_x, j_y, k_y;
    Quadrature2D<T, Quad, Integral2D<T,Quad,BasisFunc_x,BasisFunc_y > > quadrature;

	Integral2D(const Function2D<T> &_F, const BasisFunc_x &_phi_x, const BasisFunc_y &_phi_y);

    T
    operator()(int _j_x, int _k_x, int _j_y, int _k_y) const;

    T
    integrand(T x, T y) const;
};

} // namespace lawa

#include <lawa/righthandsides/integrals2d.tcc>


#endif	//LAWA_INTEGRALS2D_H
