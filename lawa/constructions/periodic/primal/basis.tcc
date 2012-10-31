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

namespace lawa {

template <typename T>
Basis<T,Primal,Periodic,CDF>::Basis(int _d, int _d_, int j)
    : d(_d), d_(_d_), j0(j), mra(d,d_,j), mra_(d,d_,j), 
      psi(d,d_), M1(psi), _j(j), refinementbasis(_d, _d_, j)
{
	if(d == 2 && d_ == 2){
        _periodicRefCoeffs = new DenseVector<Array<long double> >[1];
        _periodicRefCoeffs[0].engine().resize(5,0);
        _periodicRefCoeffs[0] = 1.L/(4.L*std::sqrt(2.L)), 1.L/(2.L*std::sqrt(2.L)), - 3.L/(2.L*std::sqrt(2.L)), 1.L/(2.L*std::sqrt(2.L)), 1.L/(4.L*std::sqrt(2.L));

        _rightRefCoeffs = new DenseVector<Array<long double> >[2];
        _rightRefCoeffs[0].engine().resize(6,0);
        _rightRefCoeffs[0] =  1.L/(4.L*std::sqrt(2.L)), 1.L/(2.L*std::sqrt(2.L)), - 3.L/(2.L*std::sqrt(2.L)),
        					  1.L/(2.L*std::sqrt(2.L)), 1.L/(2.L*std::sqrt(2.L)) + 1.L/(4.L*std::sqrt(2.L)), 1.L/(4.L*std::sqrt(2.L));

        _rightRefCoeffs[1].engine().resize(6,0);
        _rightRefCoeffs[1] =  1.L/(4.L*std::sqrt(2.L)), 1.L/(2.L*std::sqrt(2.L)), 1.L/(2.L*std::sqrt(2.L)) - 3.L/(2.L*std::sqrt(2.L)), - 3.L/(2.L*std::sqrt(2.L)), 1.L/(2.L*std::sqrt(2.L)), 1.L/(4.L*std::sqrt(2.L));

		_split = new long[2];
		_split[0] = 5;
		_split[1] = 3;

	}
}

template <typename T>
Basis<T,Primal,Periodic, CDF>::~Basis()
{
	delete[] _periodicRefCoeffs;
	/*
	 * TODO
	 */
}

template <typename T>
int
Basis<T,Primal,Periodic,CDF>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Primal,Periodic,CDF>::setLevel(int j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
const BasisFunction<T,Primal,Periodic,CDF> &
Basis<T,Primal,Periodic,CDF>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi;
    } else {
        return psi;
    }
}

template <typename T>
int
Basis<T,Primal,Periodic,CDF>::cardJ(int j) const
{
    assert(j>=j0);

    return pow2i<T>(j);
}

template <typename T>
int
Basis<T,Primal,Periodic,CDF>::cardJL(int /*j*/) const
{
	return std::max( std::ceil((d + d_ - 2)/2.0 - 1), 0.0);
}

template <typename T>
int
Basis<T,Primal,Periodic,CDF>::cardJI(int j) const
{
    assert(j>=j0);

    return std::max(pow2i<T>(j) - cardJL() - cardJR(), 0.0);
}

template <typename T>
int
Basis<T,Primal,Periodic,CDF>::cardJR(int /*j*/) const
{
    return std::ceil((d + d_)/2.0 - 1) + 1;
}

template <typename T>
const Range<int>
Basis<T,Primal,Periodic,CDF>::rangeJ(int j) const
{
    assert(j>=j0);

    return Range<int>(1,pow2i<T>(j));
}


template <typename T>
const Range<int>
Basis<T,Primal,Periodic,CDF>::rangeJL(int j) const
{
    assert(j>=j0);

    return Range<int>(1,cardJL());
}


template <typename T>
const Range<int>
Basis<T,Primal,Periodic,CDF>::rangeJI(int j) const
{
    assert(j>=j0);

    return Range<int>(cardJL() + 1,pow2i<T>(j) - cardJR());
}

template <typename T>
const Range<int>
Basis<T,Primal,Periodic,CDF>::rangeJR(int j) const
{
    assert(j>=j0);

    return Range<int>(pow2i<T>(j) - std::ceil((d + d_)/2.0 - 1), pow2i<T>(j));
}


} // namespace lawa

