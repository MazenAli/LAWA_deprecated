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

#ifndef LAWA_ADAPTIVE_REFERENCESOLUTIONS_REFERENCESOLUTIONS3D_H
#define LAWA_ADAPTIVE_REFERENCESOLUTIONS_REFERENCESOLUTIONS3D_H 1

namespace lawa {

template <typename T, typename Basis3D, typename BilinearForm>
struct ReferenceSolutionTensor3D
{
};

template<typename T, typename Basis3D>
struct ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >
{
    static int nr;
    static T c;
    static DomainType domain1, domain2, domain3;

    static DenseVector<Array<T> > sing_pts_x, sing_pts_y, sing_pts_z;	//aligned singularities

    static flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > deltas_x, deltas_y, deltas_z;

    static void
    setExample(int _nr, const HelmholtzOperator3D<T,Basis3D> &a, DomainType domain1, DomainType domain2, DomainType domain3);

    static T
    exact(T x, T y, T z);

    static T
    exact_dx(T x, T y, T z);

    static T
    exact_dy(T x, T y, T z);

    static T
    exact_dz(T x, T y, T z);

    static T
    exact_x(T x);

    static T
    exact_x(T x, int deriv_x);

    static T
    exact_y(T y);

    static T
    exact_y(T y, int deriv_y);

    static T
    exact_z(T z);

    static T
    exact_z(T z, int deriv_z);

    static T
    rhs_x(T x);

    static T
    rhs_y(T y);

    static T
    rhs_z(T z);

    static T
    H1norm();
};



}	//namespace lawa

#include <lawa/adaptive/referencesolutions/referencesolutions3d.tcc>


#endif // LAWA_ADAPTIVE_REFERENCESOLUTIONS_REFERENCESOLUTIONS3D_H