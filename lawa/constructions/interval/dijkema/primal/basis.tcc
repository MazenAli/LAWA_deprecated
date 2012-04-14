/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2011  Mario Rometsch, Alexander Stippler.

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
#include <lawa/constructions/interval/initial_stable_completion.h>

namespace lawa {

template <typename T>
Basis<T,Primal,Interval,Dijkema>::Basis(int _d, int _d_, int j)
    : mra(_d, j), mra_(_d, _d_, j),
      d(_d), d_(_d_), mu(d&1),
      min_j0(mra_.min_j0), j0(mra_.j0), _bc(2,0), _j(-1), psi(*this), refinementbasis(*this)
{
    GeMatrix<FullStorage<T,ColMajor> > Mj1, Mj1_;
    initial_stable_completion(mra,mra_,Mj1,Mj1_);
    const int cons_j = ((d==2) && ((d_==2)||(d_==4))) ? mra_.min_j0+1 : mra_.min_j0;
    M1 = RefinementMatrix<T,Interval,Dijkema>(d+d_-2, d+d_-2, Mj1, min_j0, cons_j);
    _j = std::max(min_j0,j);
    setLevel(_j);

    if (d==2 && d_==2) {
        // left part
        _leftRefCoeffs = new DenseVector<Array<long double> >[2];
        _leftRefCoeffs[0].engine().resize(3,0);
        _leftRefCoeffs[0] =  1.2374368670764579, -0.3535533905932737, -0.1767766952966368;
        _leftRefCoeffs[1].engine().resize(5,0);
        _leftRefCoeffs[1] =  -0.1767766952966368, -0.3535533905932737, 1.0606601717798211, -0.3535533905932737, -0.1767766952966368;
        _leftOffsets = new long[2];
        _leftOffsets[0] =  2;
        _leftOffsets[1] =  2;

        // inner part
        _innerRefCoeffs = new DenseVector<Array<long double> >[1];
        _innerRefCoeffs[0].engine().resize(5,0);
        _innerRefCoeffs[0] =  -0.1767766952966368, -0.3535533905932737, 1.0606601717798211, -0.3535533905932737, -0.1767766952966368;
        _innerOffsets = new long[1];
        _innerOffsets[0] =  -2;

        // inner part
        _rightRefCoeffs = new DenseVector<Array<long double> >[2];
        _rightRefCoeffs[0].engine().resize(5,0);
        _rightRefCoeffs[0] =  -0.1767766952966368, -0.3535533905932737, 1.0606601717798211, -0.3535533905932737, -0.1767766952966368;
        _rightRefCoeffs[1].engine().resize(3,0);
        _rightRefCoeffs[1] =  -0.1767766952966368, -0.3535533905932737, 1.2374368670764579;
        _rightOffsets = new long[2];
        _rightOffsets[0] =  - 2;
        _rightOffsets[1] =    0;
    }
}

template <typename T>
int
Basis<T,Primal,Interval,Dijkema>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Primal,Interval,Dijkema>::setLevel(int j) const
{
    assert(j>=min_j0);
    _j = j;
    M1.setLevel(_j);
    mra.setLevel(_j);
    mra_.setLevel(_j);
}

template <typename T>
template <BoundaryCondition BC>
void
Basis<T,Primal,Interval,Dijkema>::enforceBoundaryCondition()
{
    if ((_bc(0)==0) && (_bc(1)==0)) {
        _bc(0) = _bc(1) = 1;
        mra.enforceBoundaryCondition<BC>();
        mra_.enforceBoundaryCondition<BC>();
        GeMatrix<FullStorage<T,ColMajor> > Mj1, Mj1_;
        initial_stable_completion(mra,mra_,Mj1,Mj1_);
        const int cons_j = ((d==2) && ((d_==2)||(d_==4))) ? mra_.min_j0+1 : mra_.min_j0;
        M1 = RefinementMatrix<T,Interval,Dijkema>(d+d_-2, d+d_-2, Mj1, min_j0, cons_j);
        setLevel(_j);
    }
    // Refinement coefficients only in double prec. available due to missing support of higher
    // precision in blas routines.
    if (d==2 && d_==2) {
        // left part
        _leftRefCoeffs = new DenseVector<Array<long double> >[2];
        _leftRefCoeffs[0].engine().resize(3,0);
        _leftRefCoeffs[0] =  1.2374368670764579, -0.3535533905932737, -0.1767766952966368;
        _leftRefCoeffs[1].engine().resize(5,0);
        _leftRefCoeffs[1] =  -0.1767766952966368, -0.3535533905932737, 1.0606601717798211, -0.3535533905932737, -0.1767766952966368;
        _leftOffsets = new long[2];
        _leftOffsets[0] =  2;
        _leftOffsets[1] =  2;

        // inner part
        _innerRefCoeffs = new DenseVector<Array<long double> >[1];
        _innerRefCoeffs[0].engine().resize(5,0);
        _innerRefCoeffs[0] =  -0.1767766952966368, -0.3535533905932737, 1.0606601717798211, -0.3535533905932737, -0.1767766952966368;
        _innerOffsets = new long[1];
        _innerOffsets[0] =  -2;

        // inner part
        _rightRefCoeffs = new DenseVector<Array<long double> >[2];
        _rightRefCoeffs[0].engine().resize(5,0);
        _rightRefCoeffs[0] =  -0.1767766952966368, -0.3535533905932737, 1.0606601717798211, -0.3535533905932737, -0.1767766952966368;
        _rightRefCoeffs[1].engine().resize(3,0);
        _rightRefCoeffs[1] =  -0.1767766952966368, -0.3535533905932737, 1.2374368670764579;
        _rightOffsets = new long[2];
        _rightOffsets[0] =  - 2;
        _rightOffsets[1] =    0;
    }
    else if (d==3 && d_==3) {
        // left part
        _leftRefCoeffs = new DenseVector<Array<long double> >[4];
        _leftRefCoeffs[0].engine().resize(7,0);
        _leftRefCoeffs[0] =  0.825002583178343, -0.31011594165461, -0.063691491010539, 0.044358015602397, 0.014032578184197, -0.000847605393676, -0.0002825351312254;
        _leftRefCoeffs[1].engine().resize(7,0);
        _leftRefCoeffs[1] = -0.084857408555486,  0.12118335399876,  0.702979696218226757, -0.707687538747675, -0.110818350898946, 0.1407121822690638, 0.0469040607563546;
        _leftRefCoeffs[2].engine().resize(8,0);
        _leftRefCoeffs[2] = -0.046875, -0.140625, 0.109375, 0.703125, -0.703125, -0.109375, 0.140625, 0.046875;
        _leftRefCoeffs[3].engine().resize(8,0);
        _leftRefCoeffs[3] = -0.046875, -0.140625, 0.109375, 0.703125, -0.703125, -0.109375, 0.140625, 0.046875;
        _leftOffsets = new long[4];
        _leftOffsets[0] =  2;
        _leftOffsets[1] =  2;
        _leftOffsets[2] =  3;
        _leftOffsets[3] =  5;

        // inner part
        _innerRefCoeffs = new DenseVector<Array<long double> >[2];
        _innerRefCoeffs[0].engine().resize(8,0);
        _innerRefCoeffs[0] =  0.046875,  0.140625, -0.109375, -0.703125, 0.703125, 0.109375, -0.140625, -0.046875;
        _innerRefCoeffs[1].engine().resize(8,0);
        _innerRefCoeffs[1] = -0.046875, -0.140625, 0.109375, 0.703125, -0.703125, -0.109375, 0.140625, 0.046875;
        _innerOffsets = new long[2];
        _innerOffsets[0] =  -3;
        _innerOffsets[1] =  -3;

        // inner part
        _rightRefCoeffs = new DenseVector<Array<long double> >[4];
        _rightRefCoeffs[0].engine().resize(8,0);
        _rightRefCoeffs[0] =  0.046875, 0.140625, -0.109375, -0.703125, 0.703125, 0.109375, -0.140625, -0.046875;
        _rightRefCoeffs[1].engine().resize(8,0);
        _rightRefCoeffs[1] =  0.046875, 0.140625, -0.109375, -0.703125, 0.703125, 0.109375, -0.140625, -0.046875;
        _rightRefCoeffs[2].engine().resize(7,0);
        _rightRefCoeffs[2] = 0.0469040607563546, 0.1407121822690638,  -0.110818350898946,  -0.707687538747675, 0.702979696218226757, 0.12118335399876, -0.084857408555486;
        _rightRefCoeffs[3].engine().resize(7,0);
        _rightRefCoeffs[3] = -0.0002825351312254, -0.000847605393676, 0.014032578184197, 0.044358015602397, -0.063691491010539, -0.31011594165461, 0.825002583178343;
        _rightOffsets = new long[4];
        _rightOffsets[0] =  - 7;
        _rightOffsets[1] =  - 5;
        _rightOffsets[2] =  - 3;
        _rightOffsets[3] =  - 3;
    }
}

template <typename T>
const BasisFunction<T,Primal,Interval,Dijkema> &
Basis<T,Primal,Interval,Dijkema>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi; 
    } else {
        return psi;
    }
}

// cardinalities of whole, left, inner, right index sets (primal).
template <typename T>
int
Basis<T,Primal,Interval,Dijkema>::cardJ(int j) const
{
    assert(j>=min_j0);
    return pow2i<T>(j);
}

template <typename T>
int
Basis<T,Primal,Interval,Dijkema>::cardJL(int j) const
{
    assert(j>=min_j0);
    return d + d_ - 2;
}

template <typename T>
int
Basis<T,Primal,Interval,Dijkema>::cardJI(int j) const
{
    assert(j>=min_j0);
    return pow2i<T>(j) - 2*(d + d_ - 2);
}

template <typename T>
int
Basis<T,Primal,Interval,Dijkema>::cardJR(int j) const
{
    assert(j>=min_j0);
    return d + d_ - 2;
}

// ranges of whole, left, inner, right index sets (primal).
template <typename T>
const Range<int>
Basis<T,Primal,Interval,Dijkema>::rangeJ(int j) const
{
    assert(j>=min_j0);
    return _(1,pow2i<T>(j));
}

template <typename T>
const Range<int>
Basis<T,Primal,Interval,Dijkema>::rangeJL(int j) const
{
    assert(j>=min_j0);
    return _(1,d+d_-2);
}

template <typename T>
const Range<int>
Basis<T,Primal,Interval,Dijkema>::rangeJI(int j) const
{
    assert(j>=min_j0);
    return _(d+d_-1,pow2i<T>(j)-(d+d_-2));
}

template <typename T>
const Range<int>
Basis<T,Primal,Interval,Dijkema>::rangeJR(int j) const
{
    assert(j>=min_j0);
    return _(pow2i<T>(j)-(d+d_-3),pow2i<T>(j));
}

} // namespace lawa

