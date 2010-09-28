/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009 Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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


#ifndef LAWA_OPERATORS_PRECONDITIONER_H
#define LAWA_OPERATORS_PRECONDITIONER_H 1

#include <lawa/adaptive/index.h>
#include <lawa/operators/operators.h>

namespace lawa {
/*  A preconditioner for a 1D problem must have the following form:
 *
template <typename T, typename Index, typename Basis, typename BilinearForm>
class Preconditioner
{
public:
    Preconditioner(const BilinearForm &basis);

    T
    operator()(XType xtype1, int j1, int k1) const;

    T
    operator()(const Index &index) const;

};
*/

template <typename T>
class NoPreconditioner1D
{

public:
    NoPreconditioner1D(void) {

    };

    T
    operator()(XType xtype1, int j1, int k1) const;

    T
    operator()(const Index1D &index) const;

};

template <typename T, typename Basis, typename BilinearForm>
class DiagonalMatrixPreconditioner1D
{

public:
    DiagonalMatrixPreconditioner1D(const BilinearForm &a);

    T
    operator()(XType xtype1, int j1, int k1) const;

    T
    operator()(const Index1D &index) const;

private:
    const BilinearForm &a;
};

template <typename T, typename Basis, typename BilinearForm>
class H1Preconditioner1D
{
    typedef typename Basis::BSplineType PrimalSpline;
    typedef typename Basis::WaveletType PrimalWavelet;

public:
    H1Preconditioner1D(const BilinearForm &a);

    T
    operator()(XType xtype1, int j1, int k1) const;

    T
    operator()(const Index1D &index) const;

private:
    const Basis &basis;

    PrimalSpline phi, d_phi;
    PrimalWavelet psi, d_psi;

    Integral<T, Gauss, PrimalSpline, PrimalSpline> integral_sfsf, dd_integral_sfsf;
    Integral<T, Gauss, PrimalWavelet, PrimalWavelet> integral_ww, dd_integral_ww;

};


template <typename T, typename Basis2D, typename BilinearForm>
class DiagonalMatrixPreconditioner2D
{

public:
    DiagonalMatrixPreconditioner2D(const BilinearForm &a);

    T
    operator()(XType xtype1, int j1, int k1, XType xtype2, int j2, int k2) const;

    T
    operator()(const Index2D &index) const;

private:
    const BilinearForm &a;
};


}

#include <lawa/operators/preconditioner.tcc>


#endif /* PRECONDITIONER_H_ */
