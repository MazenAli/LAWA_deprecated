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

template <typename T>
MRA<T,Primal,Periodic,CDF>::MRA(int _d, int _d_, int j)
    : d(_d), d_(_d_), j0(j), M0(phi), _j(j), mu(d&1), phi(*this)
{

    // Refinement coefficients only in double prec. available due to missing support of higher
    // precision in blas routines.
    if (d==2) {

		// inner part
		_periodicRefCoeffs = new DenseVector<Array<long double> >[1];
		_periodicRefCoeffs[0].engine().resize(3,0);
		_periodicRefCoeffs[0] =  1.L/(2.L*std::sqrt(2.L)), 1.L/std::sqrt(2.L), 1.L/(2.L*std::sqrt(2.L));
		_innerOffsets = new long[1];
        _innerOffsets[0] =  0;

        _rightRefCoeffs = new DenseVector<Array<long double> >[1];
		_rightRefCoeffs[0].engine().resize(4,0);
		_rightRefCoeffs[0] =  1.L/(2.L*std::sqrt(2.L)), 1.L/std::sqrt(2.L), 1.L/std::sqrt(2.L), 1.L/(2.L*std::sqrt(2.L));
		_split = new long[1];
		_split[0] = 2;

		/*
		 * TODO
		 */
		/*
		_innerOffsets[0] =  -2;
		_innerL2Norms = new long double[1];
		_innerL2Norms[0] =  0.;
		_innerH1SemiNorms = new long double[1];
		_innerH1SemiNorms[0] =  0.;
		*/
	}
    else if (d==3){
		// inner part
		_periodicRefCoeffs = new DenseVector<Array<long double> >[1];
        _periodicRefCoeffs[0].engine().resize(4,0);
        _periodicRefCoeffs[0] =  1.L/(4.L*std::sqrt(2.L)), 3.L/(4.L*std::sqrt(2.L)), 3.L/(4.L*std::sqrt(2.L)), 1.L/(4.L*std::sqrt(2.L));
		_innerOffsets = new long[1];
        _innerOffsets[0] =  -3;

		/*
		 * TODO
		 */
		/*
		_innerOffsets[0] =  -2;
		_innerL2Norms = new long double[1];
		_innerL2Norms[0] =  0.;
		_innerH1SemiNorms = new long double[1];
		_innerH1SemiNorms[0] =  0.;
		*/
    }

}

template <typename T>
MRA<T,Primal,Periodic,CDF>::~MRA()
{
	delete[] _periodicRefCoeffs;
	delete[] _innerOffsets;
	/*
	 * TODO
	 */
	//delete[] _innerL2Norms;
	//delete[] _innerH1SemiNorms;
}

template <typename T>
int
MRA<T,Primal,Periodic,CDF>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Primal,Periodic,CDF>::setLevel(int j) const
{
    assert(j>=0);
    _j = j;
}

//--- cardinalities of whole, left, inner, right index sets. -------------------

template <typename T>
int
MRA<T,Primal,Periodic,CDF>::cardI(int j) const
{
    assert(j>=j0);   
    return pow2i<T>(j);
}

template <typename T>
int
MRA<T,Primal,Periodic,CDF>::cardIL(int /*j*/) const
{
    return std::max( std::ceil( (d - mu)/2.0 - 1), 0.0);
}

template <typename T>
int
MRA<T,Primal,Periodic,CDF>::cardII(int j) const
{
	assert(j >= j0);
    return pow2i<T>(j) - cardIL() - cardIR();
}

template <typename T>
int
MRA<T,Primal,Periodic,CDF>::cardIR(int /*j*/) const
{
    return std::ceil((d+mu)/2.0 - 1) + 1;
}

template <typename T>
Range<int>
MRA<T,Primal,Periodic,CDF>::rangeI(int j) const
{
    assert(j>=j0);
    return Range<int>(1,pow2i<T>(j));
}

template <typename T>
Range<int>
MRA<T,Primal,Periodic,CDF>::rangeIL(int j) const
{
    assert(j>=j0);
    // If left index set is empty, return range = {1,0}
    return Range<int>(1, cardIL());
}

template <typename T>
Range<int>
MRA<T,Primal,Periodic,CDF>::rangeII(int j) const
{
    assert(j>=j0);
    return Range<int>(1 + cardIL() , pow2i<T>(j) - cardIR());
}

template <typename T>
Range<int>
MRA<T,Primal,Periodic,CDF>::rangeIR(int j) const
{
    assert(j>=j0);
    return Range<int>(pow2i<T>(j) - std::ceil((d+mu)/2.0 - 1), pow2i<T>(j));
}


} // namespace lawa

