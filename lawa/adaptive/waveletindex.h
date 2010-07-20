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

#ifndef LAWA_ADAPTIVE_WAVELETINDEX_H
#define LAWA_ADAPTIVE_WAVELETINDEX_H 1

#include <iostream>
#include <lawa/interval/dijkema/primal/basis.h>
#include <lawa/support.h>
#include <utility>

namespace lawa {
    
template <typename T, Construction Cons>
struct WaveletIndex
{
    WaveletIndex(const Basis<T,Primal,Interval,Cons> &_basis);
    WaveletIndex(const Basis<T,Primal,Interval,Cons> &_basis, int _j);
    WaveletIndex(const Basis<T,Primal,Interval,Cons> &_basis, int _j, int _k);
    WaveletIndex(const Basis<T,Primal,Interval,Cons> &_basis, int _j, int _k, XType _xtype);

    WaveletIndex & 
    operator++();

    WaveletIndex &
    operator--();
    
    int
    vectorPosition() const;
    
    Support<T>
    supportCube() const;

    Support<T>
    support() const;

    Support<T>
    descendentCube() const;

    const Basis<T,Primal,Interval,Cons> &basis;
    XType xtype;
    int j, k;
};

template <typename T, Construction Cons>
    bool
    operator==(const WaveletIndex<T,Cons> &lhs,
               const WaveletIndex<T,Cons> &rhs);

template <typename T, Construction Cons>
    bool
    operator!=(const WaveletIndex<T,Cons> &lhs, 
               const WaveletIndex<T,Cons> &rhs);

template <typename T, Construction Cons>
    std::ostream & 
    operator<<(std::ostream &out,
               const WaveletIndex<T,Cons> &index);

//------------------------------------------------------------------------------

//represents an entry of a discrete operator
template <typename T, Construction Cons>
struct Entry : public std::pair<WaveletIndex<T,Cons>,
                                WaveletIndex<T,Cons> >
{
    Entry(const WaveletIndex<T,Cons> index1,
          const WaveletIndex<T,Cons> index2);
};

template <typename T, Construction Cons>
    std::ostream & 
    operator<<(std::ostream &out, const Entry<T,Cons> &entry);

} // namespace lawa

#include <lawa/adaptive/waveletindex.tcc>

#endif // LAWA_ADAPTIVE_WAVELETINDEX_H
