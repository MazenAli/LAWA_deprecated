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

Index1D::Index1D(void)
: j(0), k(0), xtype(XBSpline), val(0)
{
}

Index1D::~Index1D(void)
{
}

Index1D::Index1D(int _j, int _k, XType _xtype)
: j(_j), k(_k), xtype(_xtype), val(xtype)
{
	/*
	std::cout << "Index1D.val = " << val << std::endl;
	for (int i=63; i>=0; --i) {
		std::cout << ((val & (1l << i)) ? 1 : 0);
	}
	std::cout<<std::endl;
	val = (val << 16) | j;
	for (int i=63; i>=0; --i) {
		std::cout << ((val & (1l << i)) ? 1 : 0);
	}
	std::cout<<std::endl;
	val = (val << 32 | (unsigned int) k);
	for (int i=63; i>=0; --i) {
		std::cout << ((val & (1l << i)) ? 1 : 0);
	}
	std::cout<<	std::endl;
	*/
	val = (((val << 16) | (unsigned short) j) << 32) | (unsigned int) k;
}

Index1D::Index1D(const Index1D &index)
: j(index.j), k(index.k), xtype(index.xtype), val(index.val)
{
}

std::ostream& operator<<(std::ostream &s, const Index1D &_i)
{
    if (_i.xtype==XBSpline) {
        s << "scaling, (" << _i.j << " , " << _i.k << ", " << _i.val << ")";
    } else {
        s << "wavelet, (" << _i.j << " , " << _i.k << ", " << _i.val << ")";
    }
    return s;
}


Index2D::Index2D(const Index1D &_index1, const Index1D &_index2)
    : index1(_index1), index2(_index2)
{
}

Index2D::~Index2D(void)
{
}

std::ostream& operator<<(std::ostream &s, const Index2D &_i)
{
    s << "(" << _i.index1 << ", " << _i.index2 << ")";
    return s;
}


template <typename Index>
Entry<Index>::Entry(const Index &_index1, const Index &_index2)
: row_index(_index1), col_index(_index2)
{
}

template <typename Index>
std::ostream& operator<<(std::ostream &s, const Entry<Index> &entry) {
    s << "[" << entry.row_index << ", " << entry.col_index  << "]";
    return s;
}

//Lexicographic implementation
/*
template <>
struct lt<Lexicographical, Index1D>
{

    inline
    bool operator()(const Index1D &left, const Index1D &right) const
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


    inline
    bool operator()(const Entry<Index1D> &left, const Entry<Index1D> &right) const
    {
        // sort Operator row-wise
        if ( !(operator()(left.row_index,right.row_index)) && !(operator()(right.row_index,left.row_index))) {
            return operator()(left.col_index, right.col_index);
        }
        else {
            return operator()(left.row_index,right.row_index);
        }
    }

};

*/

//Bitmask implementation
template <>
struct lt<Lexicographical, Index1D>
{

	inline
	bool operator()(const Index1D &left, const Index1D &right) const
	{
		return left.val < right.val;
	}

    inline
    bool operator()(const Entry<Index1D> &left, const Entry<Index1D> &right) const
    {
        // sort Operator row-wise
    	if (left.row_index.val != right.row_index.val) return left.row_index.val < right.row_index.val;
    	else										   return left.col_index.val < right.col_index.val;
    }

};

//Lexicographic implementation
/*
template <>
struct lt<Lexicographical, Index2D >
{
    inline
    bool operator()(const Index2D &left, const Index2D &right) const
    {
        if         ((left.index1.xtype==XBSpline)&&(right.index1.xtype==XWavelet)) {
            return true;
        }
        else if ((left.index1.xtype==XWavelet)&&(right.index1.xtype==XBSpline)) {
            return false;
        }
        else if ((left.index2.xtype==XBSpline)&&(right.index2.xtype==XWavelet)) {
            return true;
        }
        else if ((left.index2.xtype==XWavelet)&&(right.index2.xtype==XBSpline)) {
            return false;
        }
        else if (left.index1.j!=right.index1.j) {
            return (left.index1.j<right.index1.j);
        }
        else if (left.index1.k!=right.index1.k) {
            return (left.index1.k<right.index1.k);
        }
        else if (left.index2.j!=right.index2.j) {
            return (left.index2.j<right.index2.j);
        }
        else {
            return (left.index2.k<right.index2.k);
        }
    }

    inline
    bool operator()(const Entry<Index2D> &left, const Entry<Index2D> &right) const
    {
    // sort Operator row-wise
        if ( !(operator()(left.row_index,right.row_index)) && !(operator()(right.row_index,left.row_index))) {
            return operator()(left.col_index, right.col_index);
        }
        else {
            return operator()(left.row_index,right.row_index);
        }
    }
};
*/

template <>
struct lt<Lexicographical, Index2D >
{

	inline
	bool operator()(const Index2D &left, const Index2D &right) const
	{
		if (left.index1.val != right.index1.val) return left.index1.val < right.index1.val;
		else									 return left.index2.val < right.index2.val;
	}

    inline
    bool operator()(const Entry<Index2D> &left, const Entry<Index2D> &right) const
    {
    // sort Operator row-wise
        if (left.row_index.index1.val != right.row_index.index1.val) {
        	return left.row_index.index1.val < right.row_index.index1.val;
        }
        else if (left.row_index.index2.val != right.row_index.index2.val) {
        	return left.row_index.index2.val < right.row_index.index2.val;
        }
        if (left.col_index.index1.val != right.col_index.index1.val) {
            return left.col_index.index1.val < right.col_index.index1.val;
        }
        else {
            return left.col_index.index2.val < right.col_index.index2.val;
        }
    }

};

template <typename SortingType>
struct lt<AbsoluteValue, SortingType>
{
    bool operator()(const SortingType &left, const SortingType &right) const
        {
            return (fabs(left) > fabs(right));    //todo: Is this the right call for fabs (template?)
        }
};


} //namespace lawa
