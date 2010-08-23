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

#ifndef ADAPTIVE_AUX_PLOTTING_H
#define ADAPTIVE_AUX_PLOTTING_H 1

#include <iostream>
#include <fstream>
#include <lawa/adaptive/index.h>
#include <lawa/adaptive/indexset.h>
#include <lawa/adaptive/coefficients.h>

namespace lawa {

template <typename T, typename Basis>
void
getSingularPoints(const Basis &basis, const Coefficients<Lexicographical,T,Index1D> coeff, DenseVector<Array<T> > &sing_pts);

//Plot solution only on singular  points of the solution
template <typename T, typename Basis, typename Preconditioner>
void
plot(const Basis &basis, const Coefficients<Lexicographical,T,Index1D> coeff,
	 const Preconditioner &P, T (*u)(T), const char* filename);

//Plot solution on a fixed grid
template <typename T, typename Basis, typename Preconditioner>
void
plot(const Basis &basis, const Coefficients<Lexicographical,T,Index1D> coeff,
	 const Preconditioner &P, T (*u)(T), T a, T b, T h, const char* filename);

template <typename T, typename Basis2D, typename Preconditioner>
void
plot2D(const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D> coeff,
	   const Preconditioner &P, T (*u)(T,T), T a1, T b1, T a2, T b2, T h, const char* filename);

template <typename T, DomainType Domain, Construction Cons>
void
plotCoeff(const Coefficients<AbsoluteValue,T,Index1D > &coeff, const Basis<T,Primal,Domain,Cons> &basis, const char* filename);


}  // namespace lawa

#include <lawa/adaptive/aux/plotting.tcc>

#endif // ADAPTIVE_AUX_PLOTTING_H
