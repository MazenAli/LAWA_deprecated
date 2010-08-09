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

namespace lawa {

Index::Index(int _d, int _d_)
: d(_d), d_(_d_), j(0), k(0), xtype(XBSpline)
{

}

Index::Index(int _d, int _d_, int _j, int _k, XType _xtype)
: d(_d), d_(_d_), j(_j), k(_k), xtype(_xtype)
{

}

Index&
Index::operator= (const Index &_index)
{
    j     = _index.j;
    k     = _index.k;
    xtype = _index.xtype;
    return *this;
}

std::ostream& operator<<(std::ostream &s, const Index &_i)
{
    if (_i.xtype==XBSpline) {
        s << "scaling, (" << _i.j << " , " << _i.k << ")";
    } else {
        s << "wavelet, (" << _i.j << " , " << _i.k << ")";
    }
    return s;
}

ReallineIndex::ReallineIndex(int d, int d_)
: Index(d, d_)
{

}

ReallineIndex::ReallineIndex(int d, int d_, int j, int k, XType xtype)
: Index(d, d_, j, k, xtype)
{

}

template <typename Basis>
IntervalIndex<Basis>::IntervalIndex(const Basis &_basis, int d, int d_)
: Index(d, d_), basis(_basis)
{

}

template <typename Basis>
IntervalIndex<Basis>::IntervalIndex(const Basis &_basis, int d, int d_, int j, int k, XType xtype)
: Index(d, d_, j, k, xtype), basis(_basis)
{

}

template <typename Basis>
PeriodicIndex<Basis>::PeriodicIndex(const Basis &_basis, int d, int d_)
: Index(d, d_), basis(_basis)
{

}

template <typename Basis>
PeriodicIndex<Basis>::PeriodicIndex(const Basis &_basis, int d, int d_, int j, int k, XType xtype)
: Index(d, d_, j, k, xtype), basis(_basis)
{

}


} //namespace lawa
