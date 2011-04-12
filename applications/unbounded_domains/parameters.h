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

#ifndef APPLICATIONS_UNBOUNDEDDOMAINS_PARAMETERS_H
#define APPLICATIONS_UNBOUNDEDDOMAINS_PARAMETERS_H 1

#include <lawa/constructions/basis.h>
#include <lawa/operators/operators.h>



namespace lawa {


template <typename T, typename Basis1D, typename BilinearForm, typename Preconditioner>
class Parameters;

template <typename T>
class Parameters<T, Basis<T,Primal,R,CDF>, HelmholtzOperator1D<T,Basis<T,Primal,R,CDF> >,
                 H1Preconditioner1D<T, Basis<T,Primal,R,CDF>, HelmholtzOperator1D<T,Basis<T,Primal,R,CDF> > > >
{
public:
	const Basis<T,Primal,R,CDF> &basis;
	const HelmholtzOperator1D<T,Basis<T,Primal,R,CDF> > &Bil;
	bool w_XBSpline;
	int j0;
	T cA, CA, kappa;
	T alpha, omega, gamma, theta;	//GHSADWAV parameters (from [GHS:2007])

	Parameters(const Basis<T,Primal,R,CDF> &_basis,
			   const HelmholtzOperator1D<T,Basis<T,Primal,R,CDF> > &_Bil,
			   bool _w_XBSpline=true, int _j0=0);

	void
	getGHSADWAVParameters(T &_alpha, T &_omega, T &_gamma, T &_theta) const;

	void
	getRHS_W_XBSplineParameters(T &_left_bound, T &_right_bound, int &_J_plus_smooth,
								int &_J_plus_singular, bool &_singular_integral,
								int example) const;

	void
	getRHS_WO_XBSplineParameters(T &_left_bound, T &_right_bound, int &_J_plus_smooth, int &_J_minus_smooth,
								 int &_J_plus_singular, int &_J_minus_singular, bool &_singular_integral,
								 int example) const;


};

}	//namespace lawa

#include<applications/unbounded_domains/parameters.tcc>

#endif	//APPLICATIONS_UNBOUNDEDDOMAINS_PARAMETERS_H
