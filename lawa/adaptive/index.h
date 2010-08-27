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

#ifndef LAWA_ADAPTIVE_INDEX_H
#define LAWA_ADAPTIVE_INDEX_H 1

#include <lawa/enum.h>
#include <iostream>

namespace lawa {

class Index1D
{
public:
    int j, k;
    XType xtype;

    Index1D(void);
    ~Index1D();
    Index1D(int j, int k, XType _xtype);
    Index1D(const Index1D &index);
};


std::ostream& operator<<(std::ostream &s, const Index1D &_Index);

struct Index2D
{
    Index2D(const Index1D &index1, const Index1D &index2);
    ~Index2D();
    Index1D index1, index2;

};


std::ostream& operator<<(std::ostream &s, const Index2D &_Index);

template <typename Index>
class Entry
{
public:
    Entry(const Index &row_index, const Index &col_index);
    const Index row_index, col_index;    //todo: no copy, but only a reference possible ?!
};

template <typename Index>
std::ostream& operator<<(std::ostream &s, const Entry<Index> &entry);

template <SortingCriterion S, typename SortingType>
struct lt
{
};

} //namespace lawa

#include <lawa/adaptive/index.tcc>

#endif //LAWA_ADAPTIVE_INDEX_H
