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

#ifndef LAWA_REALLINE_PRIMAL_BASIS_H
#define LAWA_REALLINE_PRIMAL_BASIS_H 1

#include <lawa/basis.h>
#include <lawa/bspline.h>
#include <lawa/enum.h>
#include <lawa/wavelet.h>

namespace lawa {

template <typename T>
class Basis<T,Primal,R,CDF>
{
    public:
        typedef BSpline<T,Primal,R,CDF> BSplineType;
        typedef Wavelet<T,Primal,R,CDF> WaveletType;

        Basis(int _d, int _d_, int j=0);

        int 
        level() const;

        void
        setLevel(int j) const;

        const int d, d_, j0;
        MRA<T,Primal,R,CDF> mra;
        MRA<T,Dual,R,CDF> mra_;
        Wavelet<T,Primal,R,CDF> psi;
        RefinementMatrix<T,R,CDF> M1;
        
    private:
        mutable int _j;
};

} // namespace lawa

#include <lawa/realline/primal/basis.tcc>

#endif // LAWA_REALLINE_PRIMAL_BASIS_H