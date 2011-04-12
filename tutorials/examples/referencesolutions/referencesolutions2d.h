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

#ifndef TUTORIALS_EXAMPLES_REFERENCESOLUTIONS_REFERENCESOLUTIONS2D_H
#define TUTORIALS_EXAMPLES_REFERENCESOLUTIONS_REFERENCESOLUTIONS2D_H 1

namespace lawa {

template <typename T, typename Basis, typename BilinearForm>
struct ReferenceSolutionTensor2D
{
};

template<typename T, typename Basis2D>
struct ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >
{
    static int nr;
    static T c;
    static DomainType domain1, domain2;

    static DenseVector<Array<T> > sing_pts_x, sing_pts_y;	//aligned singularities

    static flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > deltas_x, deltas_y;

    static void
    setExample(int _nr, const HelmholtzOperator2D<T,Basis2D> &a, DomainType domain1, DomainType domain2);

    static T
    exact(T x, T y);

    static T
    exact_dx(T x, T y);

    static T
    exact_dy(T x, T y);

    static T
    exact_x(T x);

    static T
    exact_x(T x, int deriv_x);

    static T
    exact_y(T y);

    static T
    exact_y(T y, int deriv_y);

    static T
    rhs_x(T x);

    static T
    rhs_y(T y);

    static T
    H1norm();
};


template <typename T, typename Basis, typename BilinearForm>
struct ReferenceSolution2D
{
};

template<typename T, typename Basis2D>
struct ReferenceSolution2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >
{
    static int nr;
    static T c;
    static DomainType domain1, domain2;

    static DenseVector<Array<T> > sing_pts_x, sing_pts_y;

    static void
    setExample(int _nr, const HelmholtzOperator2D<T,Basis2D> &a, DomainType domain1, DomainType domain2);

    static T
    exact(T x, T y);

    static T
    minus_exact(T x, T y);

    static T
    exact_dx(T x, T y);

    static T
    exact_dy(T x, T y);

    static T
    rhs(T x, T y);

    static T
    exact(T x, T y, int deriv_x, int deriv_y);

    static T
    H1norm();

};


template<typename T, typename Basis2D>
struct ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >
{
	static int nr;
	static T c;
	static DomainType domain1, domain2;

	static DenseVector<Array<T> > sing_pts_t, sing_pts_x;

	static void
	setExample(int _nr, const SpaceTimeHeatOperator1D<T,Basis2D> &a, DomainType domain1, DomainType domain2);

	static T
	exact(T t, T x);
		
    static T
    dx_exact(T t, T x);

	static T
	exact_t(T t);
	
	static T
	_exact_t(T t);

	static T
	exact_t(T t, int deriv_t);

	static T
	exact_x(T x);

	static T
	exact_x(T x, int deriv_x);

    static T
    dd_exact_x(T x);

	static T
	rhs_t(T t);

	static T
	rhs_x(T x);
	
	static T
	H1_t_norm(T t);

};

}	//namespace lawa

#include <tutorials/examples/referencesolutions/referencesolutions2d.tcc>


#endif // TUTORIALS_EXAMPLES_REFERENCESOLUTIONS_REFERENCESOLUTIONS2D_H
