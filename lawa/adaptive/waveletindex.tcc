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

#include <cassert>

namespace lawa {

template <typename T, Construction Cons>
WaveletIndex<T,Cons>::WaveletIndex(const Basis<T,Primal,Interval,Cons> &_basis)
    : basis(_basis), j(basis.j0), k(basis.mra.rangeI(j).firstIndex()),
      xtype(XBSpline)
{
}

template <typename T, Construction Cons>
WaveletIndex<T,Cons>::WaveletIndex(const Basis<T,Primal,Interval,Cons> &_basis, int _j)
    : basis(_basis), j(_j), k(basis.mra.rangeI(j).firstIndex()), xtype(XBSpline)
{
    assert(j>=basis.j0);
}

template <typename T, Construction Cons>
WaveletIndex<T,Cons>::WaveletIndex(const Basis<T,Primal,Interval,Cons> &_basis, int _j, int _k)
    : basis(_basis), j(_j), k(_k), xtype(XBSpline)
{
    assert(j>=basis.j0);
    assert(k>=basis.mra.rangeI(j).firstIndex());
    assert(k<=basis.mra.rangeI(j).lastIndex());
}

template <typename T, Construction Cons>
WaveletIndex<T,Cons>::WaveletIndex(const Basis<T,Primal,Interval,Cons> &_basis, int _j , int _k,
                                        XType _xtype)
    : basis(_basis), j(_j), k(_k), xtype(_xtype)
{
    assert(j>=basis.j0);
    assert(!( (xtype==XBSpline)&&(k<basis.mra.rangeI(j).firstIndex()) ));
    assert(!( (xtype==XBSpline)&&(k>basis.mra.rangeI(j).lastIndex()) ));
    assert(!( (xtype==XWavelet)&&(k<basis.rangeJ(j).firstIndex()) ));
    assert(!( (xtype==XWavelet)&&(k>basis.rangeJ(j).lastIndex()) ));
}

template <typename T, Construction Cons>
WaveletIndex<T,Cons> &
WaveletIndex<T,Cons>::operator++()
{
    if (xtype==XBSpline) {
        if (k<basis.mra.rangeI(j).lastIndex()) {
            ++k;
        } else {
            k = basis.rangeJ(j).firstIndex();
            xtype = XWavelet;
        }
    } else {
        if (k<basis.rangeJ(j).lastIndex()) {
            ++k;
        } else {
            ++j;
            k = basis.rangeJ(j).firstIndex();
        }
    }
    return *this;
}

template <typename T, Construction Cons>
WaveletIndex<T,Cons>&
WaveletIndex<T,Cons>::operator--()
{
    if (xtype==XBSpline) {
        if (k>basis.mra.rangeI(j).firstIndex()) {
            --k;
        } else {
            assert(0);
        }
    } else {
        if (k>basis.rangeJ(j).firstIndex()) {
            --k;
        } else {
            if (j == basis.j0) {
                k = basis.mra.rangeI(j).lastIndex();
                xtype = XBSpline;
            } else {
                --j;
                k = basis.rangeJ(j).lastIndex();
            }
        }
    }
    return *this;
}

template <typename T, Construction Cons>
int
WaveletIndex<T,Cons>::vectorPosition() const
{
    if (xtype==XBSpline) {
        return k;
    } else {
        return basis.mra.rangeI(j).lastIndex() + k;
    }
}

template <typename T, Construction Cons>
Support<T>
WaveletIndex<T,Cons>::supportCube() const
{
    // (2.13) BU 2008 ACHA
    Support<T> ret;
    if (xtype==XBSpline) {
        ret.l1 = 1.0 / basis.mra.cardI(j) * (k     - basis.mra.rangeI(j).firstIndex());
        ret.l2 = 1.0 / basis.mra.cardI(j) * (k + 1 - basis.mra.rangeI(j).firstIndex());
    } else {
        ret.l1 = 1.0 / basis.cardJ(j) * (k     - basis.rangeJ(j).firstIndex());
        ret.l2 = 1.0 / basis.cardJ(j) * (k + 1 - basis.rangeJ(j).firstIndex());
    }
    return ret;
}

template <typename T, Construction Cons>
Support<T>
WaveletIndex<T,Cons>::descendentCube() const
{
    assert(Cons==DKU);
    Support<T> ret;
    if (basis.d==2 && basis.d_==2) {
        if (xtype==XBSpline) {
            int kMid=(basis.mra.cardI(j)+1)/2, kFirst=basis.mra.rangeI(j).firstIndex();
            if (k<kMid) {
                ret.l1 = 1.0 / (basis.cardI(j)+1) * (k         - kFirst);
                ret.l2 = 1.0 / (basis.cardI(j)+1) * (k + 1     - kFirst);
            } else if (k==kMid) {
                ret.l1 = 1.0 / (basis.cardI(j)+1) * (k         - kFirst);
                ret.l2 = 1.0 / (basis.cardI(j)+1) * (k + 2     - kFirst);
            } else {
                ret.l1 = 1.0 / (basis.cardI(j)+1) * (k + 1     - kFirst);
                ret.l2 = 1.0 / (basis.cardI(j)+1) * (k + 1 + 1 - kFirst);
            }
        } else {
            ret.l1 = 1.0 / basis.cardJ(j) * (k     - basis.rangeJ(j).firstIndex());
            ret.l2 = 1.0 / basis.cardJ(j) * (k + 1 - basis.rangeJ(j).firstIndex());
        }
    } else if (basis.d==3 && basis.d_==3) {
        if (xtype==XBSpline) {
            int kInnerFirst=basis.mra.cardI(j)/2,
                 kInnerLast=basis.mra.cardI(j)/2 + 3,  
                kFirst=basis.mra.rangeI(j).firstIndex();
            if (k<kInnerFirst) {
                ret.l1 = 1.0 / basis.cardJ(j) * (k         - kFirst);
                ret.l2 = 1.0 / basis.cardJ(j) * (k + 1     - kFirst);
            } else if (k>=kInnerFirst && k<=kInnerLast) {
                ret.l1 = 1.0 / basis.cardJ(j) * (k + (k-kInnerFirst-1) + 1 - kFirst);
                ret.l2 = 1.0 / basis.cardJ(j) * (k + (k-kInnerFirst-1) + 3 - kFirst);
            } else {
                ret.l1 = 1.0 / basis.cardJ(j) * (k + 4     - kFirst);
                ret.l2 = 1.0 / basis.cardJ(j) * (k + 5     - kFirst);
            }
        } else {
            ret.l1 = 1.0 / basis.cardJ(j) * (k     - basis.rangeJ(j).firstIndex());
            ret.l2 = 1.0 / basis.cardJ(j) * (k + 1 - basis.rangeJ(j).firstIndex());
        }
    } else {
        std::cerr << "Implement descendentCube() function for d=" << basis.d
            << ", d_=" << basis.d_ << " first!" << std::endl;
        assert(false);
    }
    return ret;
}

template <typename T, Construction Cons>
Support<T>
WaveletIndex<T,Cons>::support() const
{
    if (xtype==XBSpline) {
        return basis.mra.phi.support(j,k);
    } else {
        return basis.psi.support(j,k);
    }
}

template <typename T, Construction Cons>
std::ostream &
operator<<(std::ostream &out, const WaveletIndex<T,Cons> &index)
{
    if (index.xtype==XBSpline) {
        out << "scaling, (" << index.j << " , " << index.k << ")";
    } else {
        out << "wavelet, (" << index.j << " , " << index.k << ")";
    }
    return out;
}

template <typename T, Construction Cons>
bool
operator==(const WaveletIndex<T,Cons> &lhs,
           const WaveletIndex<T,Cons> &rhs)
{
    return (lhs.j==rhs.j) && (lhs.k==rhs.k) && (lhs.xtype==rhs.xtype);
}

template <typename T, Construction Cons>
bool
operator!=(const WaveletIndex<T,Cons> &lhs, 
           const WaveletIndex<T,Cons> &rhs)
{
    return (lhs.j!=rhs.j || lhs.k!=rhs.k || lhs.xtype!=rhs.xtype);
}

template <typename T, Construction Cons>
Entry<T,Cons>::Entry(const WaveletIndex<T,Cons> first,
                     const WaveletIndex<T,Cons> second)
    : std::pair<WaveletIndex<T,Cons>,WaveletIndex<T,Cons> >(first,second)
{
}

template <typename T, Construction Cons>
std::ostream& operator<< (std::ostream &out, const Entry<T,Cons> &entry)
{
    out << "entry: (" << entry.first << ", " << entry.second << ")";
    return out;
}

} // namespace lawa
