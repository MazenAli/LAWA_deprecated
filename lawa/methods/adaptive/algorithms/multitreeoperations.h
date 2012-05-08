/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Kristina Steih, Alexander Stippler.

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

#ifndef  LAWA_METHODS_ADAPTIVE_ALGORITHMS_MULTITREEOPERATIONS_H
#define  LAWA_METHODS_ADAPTIVE_ALGORITHMS_MULTITREEOPERATIONS_H 1


#include <iostream>
#include <lawa/constructions/constructions.h>
#include <lawa/methods/adaptive/datastructures/datastructures.h>

namespace lawa {

template <typename T, typename Basis>
void
extendMultiTree(const Basis &basis, const Coefficients<Lexicographical,T,Index2D>  &v,
                Coefficients<Lexicographical,T,Index2D>  &C_v);

template <typename T, typename Basis>
void
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v);

template <typename T, typename Basis>
void
getSparseGridVector(const Basis &basis, Coefficients<Lexicographical,T,Index2D> &v, int j, T gamma);

/*
template <typename Index, typename Basis>
void
extendMultiTree(const Basis &basis, const Index &index2d, IndexSet<Index> &Lambda);

template <typename Index, typename Basis>
void
extendMultiTree2(const Basis &basis, const Index &index2d, const int offset, IndexSet<Index> &Lambda);
*/

}

#include <lawa/methods/adaptive/algorithms/multitreeoperations.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_MULTITREEOPERATIONS_H
