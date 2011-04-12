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

#ifndef LAWA_RIGHTHANDSIDES_SEPARABLERHS3D_H
#define LAWA_RIGHTHANDSIDES_SEPARABLERHS3D_H 1

#include <lawa/functiontypes/separablefunctionnd.h>
#include <lawa/integrals/integrals.h>

namespace lawa {

template<typename T, typename Basis3D>
class SeparableRHS3D
{
    private:
        const Basis3D& basis;
        const SeparableFunction3D<T>& F;
        const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_x;
	    const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_y;
	    const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_z;

        typedef typename Basis3D::FirstBasisType::BSplineType PrimalSpline_x;
        typedef typename Basis3D::SecondBasisType::BSplineType PrimalSpline_y;
        typedef typename Basis3D::ThirdBasisType::BSplineType PrimalSpline_z;
        typedef typename Basis3D::FirstBasisType::WaveletType PrimalWavelet_x;
        typedef typename Basis3D::SecondBasisType::WaveletType PrimalWavelet_y;
        typedef typename Basis3D::ThirdBasisType::WaveletType PrimalWavelet_z;

        PrimalSpline_x phi_x;
        PrimalSpline_y phi_y;
        PrimalSpline_y phi_z;
        PrimalWavelet_x psi_x;
        PrimalWavelet_y psi_y;
        PrimalWavelet_y psi_z;

        Integral<T, Gauss, PrimalSpline_x, Function<T> > integral_sff_x;
        Integral<T, Gauss, PrimalSpline_y, Function<T> > integral_sff_y;
        Integral<T, Gauss, PrimalSpline_y, Function<T> > integral_sff_z;
        Integral<T, Gauss, PrimalWavelet_x, Function<T> > integral_wf_x;
        Integral<T, Gauss, PrimalWavelet_y, Function<T> > integral_wf_y;
        Integral<T, Gauss, PrimalWavelet_y, Function<T> > integral_wf_z;

    public:
        SeparableRHS3D(const Basis3D& _basis, const SeparableFunction3D<T>& _F,
        	const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_x,
        	const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_y,
        	const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas_z, int order);

        T
        operator()(XType XisSpline, int j_x, int k_x,
                   XType YisSpline, int j_y, int k_y,
                   XType ZisSpline, int j_z, int k_z) const;

        T
        operator()(const Index3D &index) const;
};

} // namespace lawa

#include <lawa/righthandsides/separablerhs3d.tcc>

#endif // LAWA_RIGHTHANDSIDES_SEPARABLERHS3D_H
