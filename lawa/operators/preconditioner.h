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

template <typename T, typename Index>
class NoPreconditioner
{

public:
    NoPreconditioner(void) {

    };

    T
    operator()(const Index &index) const;

};

template <typename T>
class DiagonalLevelPreconditioner1D
{
public:
    DiagonalLevelPreconditioner1D(void) { };

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


/*
 * 2d preconditioner
 */

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

template <typename T, typename Basis2D, typename BilinearForm>
class RightNormPreconditioner2D
{
};

template <typename T, typename Basis2D>
class RightNormPreconditioner2D<T,Basis2D, SpaceTimeHeatOperator1D<T, Basis2D> >
{
private:
    const SpaceTimeHeatOperator1D<T, Basis2D> &a;

    typedef typename Basis2D::FirstBasisType::BSplineType PrimalSpline_t;
    typedef typename Basis2D::SecondBasisType::BSplineType PrimalSpline_x;
    typedef typename Basis2D::FirstBasisType::WaveletType PrimalWavelet_t;
    typedef typename Basis2D::SecondBasisType::WaveletType PrimalWavelet_x;

    PrimalSpline_t phi_t, d_phi_t;
    PrimalSpline_x phi_x, d_phi_x;
    PrimalWavelet_t psi_t, d_psi_t;
    PrimalWavelet_x psi_x, d_psi_x;

    Integral<T, Gauss, PrimalSpline_t, PrimalSpline_t>      integral_sfsf_t,
                                                         dd_integral_sfsf_t;
    Integral<T, Gauss, PrimalSpline_x, PrimalSpline_x>      integral_sfsf_x,
                                                         dd_integral_sfsf_x;
    Integral<T, Gauss, PrimalWavelet_t, PrimalWavelet_t>    integral_ww_t,
                                                         dd_integral_ww_t;
    Integral<T, Gauss, PrimalWavelet_x, PrimalWavelet_x>    integral_ww_x,
                                                         dd_integral_ww_x;

public:
    RightNormPreconditioner2D(const SpaceTimeHeatOperator1D<T, Basis2D> &a);

    T
    operator()(XType xtype1, int j1, int k1, XType xtype2, int j2, int k2) const;

    T
    operator()(const Index2D &index) const;
};

template <typename T, typename Basis2D, typename BilinearForm>
class LeftNormPreconditioner2D
{

};

template <typename T, typename Basis2D>
class LeftNormPreconditioner2D<T,Basis2D, SpaceTimeHeatOperator1D<T, Basis2D> >
{
private:
    const SpaceTimeHeatOperator1D<T, Basis2D> &a;

    typedef typename Basis2D::SecondBasisType::BSplineType PrimalSpline_x;
    typedef typename Basis2D::SecondBasisType::WaveletType PrimalWavelet_x;

    PrimalSpline_x  phi_x, d_phi_x;
    PrimalWavelet_x psi_x, d_psi_x;

    Integral<T, Gauss, PrimalSpline_x, PrimalSpline_x>      integral_sfsf_x,
                                                         dd_integral_sfsf_x;
    Integral<T, Gauss, PrimalWavelet_x, PrimalWavelet_x>    integral_ww_x,
                                                         dd_integral_ww_x;

public:
    LeftNormPreconditioner2D(const SpaceTimeHeatOperator1D<T, Basis2D> &a);

    T
    operator()(XType xtype1, int j1, int k1, XType xtype2, int j2, int k2) const;

    T
    operator()(const Index2D &index) const;
};



/*
 * 3d preconditioner
 */

template <typename T, typename Basis3D, typename BilinearForm>
class DiagonalMatrixPreconditioner3D
{

public:
    DiagonalMatrixPreconditioner3D(const BilinearForm &a);

    T
    operator()(XType xtype1, int j1, int k1, XType xtype2, int j2, int k2,
    		   XType xtype3, int j3, int k3) const;

    T
    operator()(const Index3D &index) const;

private:
    const BilinearForm &a;
};


}

#include <lawa/operators/preconditioner.tcc>


#endif /* PRECONDITIONER_H_ */
