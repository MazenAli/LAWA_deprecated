/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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

#ifndef APPLICATIONS_CANUTOPROJECT_LINEARTENSORINTERPOLATIONPIC2D_H
#define APPLICATIONS_CANUTOPROJECT_LINEARTENSORINTERPOLATIONPIC2D_H 1

#include <lawa/settings/enum.h>
#include <lawa/flensforlawa.h>

namespace lawa {

template<typename T>
struct LinearTensorInterpolationPic2D
{
    typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
    typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
    typedef Basis<T,Primal,R,CDF>                                       PrimalBasis;

    static int N1, N2;
    static DenseVectorT sing_pts_x, sing_pts_y;
    static DenseMatrixT coeffs;

    static PrimalBasis  linearTensorInterpolBasis;

    static void
    readPicture(const char* filename);

    static T
    evaluateInterpolation(T x1, T x2);

    static void
    plotInterpolation(const char* filename, T h1, T h2);


};


} // namespace lawa

#include <applications/canutoproject/lineartensorinterpolationpic2d.tcc>


#endif // APPLICATIONS_CANUTOPROJECT_LINEARTENSORINTERPOLATIONPIC2D_H
