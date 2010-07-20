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

#ifndef LAWA_ADAPTIVE_EVALUATE_H
#define LAWA_ADAPTIVE_EVALUATE_H 1

#include <lawa/enum.h>
#include <lawa/adaptive/waveletcoefficient.h>

namespace lawa {

//--- evaluate adaptive wavelet coefficients --------------

template <typename T,Construction Cons>
    T
    evaluate(const Coefficient<Uniform,T,Cons> &coeff, T x,
             int deriv = 0);

template <typename T,Construction Cons>
    T
    evaluate(const Coefficient<Lexicographical,T,Cons> &coeff, T x,
             int deriv = 0);

template <typename T,Construction Cons>
    DenseVector<Array<T> >
    evaluate(const Coefficient<Lexicographical,T,Cons> &coeff, 
             const DenseVector<Array<T> > &x, int deriv = 0);

template <typename T,Construction Cons>
    T
    evaluate(const Coefficient<AbsoluteValue,T,Cons> &coeff, T x,
             int deriv = 0);

} // namespace lawa

#include <lawa/adaptive/evaluate.tcc>

#endif // LAWA_ADAPTIVE_EVALUATE_H