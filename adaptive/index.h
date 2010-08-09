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

#ifndef INDEX_H
#define INDEX_H 1

#include <lawa/enum.h>

namespace lawa {

class Index
{
public:
	const int d, d_;
	int j, k;
	XType xtype;

	Index(const int d, const int d_);
    Index(const int d, const int d_, int j, int k, XType _xtype);

    Index&
    operator= (const Index &_index);

};

std::ostream& operator<<(std::ostream &s, const Index &_Index);

class ReallineIndex : public Index
{
public:
	ReallineIndex(int d, int d_);
	ReallineIndex(int d, int d_, int j, int k, XType _xtype);
};


template <typename Basis>
class IntervalIndex : public Index
{
public:
	const Basis &basis;

	IntervalIndex(const Basis &_basis, int d, int d_);
	IntervalIndex(const Basis &_basis, int d, int d_, int j, int k, XType _xtype);
};

template <typename Basis>
class PeriodicIndex : public Index
{
public:
	const Basis &basis;

	PeriodicIndex(const Basis &_basis, int d, int d_);
	PeriodicIndex(const Basis &_basis, int d, int d_, int j, int k, XType _xtype);
};


} //namespace lawa

#include "index.tcc"

#endif
