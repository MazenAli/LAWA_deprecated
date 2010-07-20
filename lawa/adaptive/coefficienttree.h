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

#ifndef LAWA_ADAPTIVE_COEFFICIENTTREE_H
#define LAWA_ADAPTIVE_COEFFICIENTTREE_H

#include <lawa/adaptive/tree.h>

namespace lawa {

template <typename T, Construction Cons >
bool
isDescendant(const WaveletIndex<T,Cons> &i1,
             const WaveletIndex<T,Cons> &i2);

template <typename T,Construction Cons >
Coefficient<Lexicographical,T,Cons>
completeTree(const Coefficient<Lexicographical,T,Cons> &coeff);

template <typename T,Construction Cons >
tree<WaveletIndex<T,Cons> >
convert(const Coefficient<Lexicographical,T,Cons> &coeff);

template <typename T,Construction Cons >
void
addToTree(tree<WaveletIndex<T,Cons> > &tr, const WaveletIndex<T,Cons> &index);

template <typename T,Construction Cons >
typename tree<WaveletIndex<T,Cons> >::iterator
findInTree(const tree<WaveletIndex<T,Cons> > &tr, const WaveletIndex<T,Cons> &index);

template <typename T,Construction Cons>
std::ostream&
operator<<(std::ostream &s, const tree<WaveletIndex<T,Cons> > &_tr);

} // namespace lawa

#include <lawa/adaptive/coefficienttree.tcc>

#endif
