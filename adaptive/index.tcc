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

Index1d::Index1d(void)
: j(0), k(0), xtype(XBSpline)
{

}

Index1d::Index1d(int _j, int _k, XType _xtype)
: j(_j), k(_k), xtype(_xtype)
{

}

std::ostream& operator<<(std::ostream &s, const Index1d &_i)
{
    if (_i.xtype==XBSpline) {
        s << "scaling, (" << _i.j << " , " << _i.k << ")";
    } else {
        s << "wavelet, (" << _i.j << " , " << _i.k << ")";
    }
    return s;
}


template <>
struct lt<Lexicographical, Index1d>
{
    inline
    bool operator()(const Index1d &left, const Index1d &right) const
    {
        if (left.j!=right.j) {
            return (left.j<right.j);
        } else {
            if ((left.xtype==XBSpline)&&(right.xtype==XWavelet)) {
                return true;
            } else if ((left.xtype==XWavelet)&&(right.xtype==XBSpline)) {
                return false;
            } else {
                return (left.k<right.k);
            }
        }
    }
};

template <typename SortingType>
struct lt<AbsoluteValue, SortingType>
{
    bool operator()(const SortingType &left, const SortingType &right) const
        {
            return (fabs(left) > fabs(right));	//todo: Is this the right call for fabs (template?)
        }
};


} //namespace lawa
