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

namespace lawa {

template <typename T, QuadratureType Quad, typename BasisFunc_x, typename BasisFunc_y>
T
_integrate(const Integral2D<T,Quad,BasisFunc_x,BasisFunc_y> integral) {
	const BasisFunc_x &phi_x   = integral.phi_x;
	const BasisFunc_y &phi_y   = integral.phi_y;

	// merge singular points of bspline/wavelet and function to one list.
	// -> implicit assumption: second.singularPoints are sorted!!
	DenseVector<Array<T> > SingularPoints_phi_x
	                              = phi_x.singularSupport(integral.j_x,integral.k_x);
	int n_phi_x = SingularPoints_phi_x.length(),
	    n_F_x   = integral.F.singularPts_x.length();

	DenseVector<Array<T> > AllSingularPoints_x(n_phi_x + n_F_x);

	std::merge(SingularPoints_phi_x.engine().data(),
			   SingularPoints_phi_x.engine().data() + n_phi_x,
	           integral.F.singularPts_x.engine().data(),
	           integral.F.singularPts_x.engine().data() + n_F_x,
	           AllSingularPoints_x.engine().data());

	DenseVector<Array<T> > SingularPoints_phi_y
		                               = phi_y.singularSupport(integral.j_y,integral.k_y);
	int n_phi_y = SingularPoints_phi_y.length(),
		n_F_y   = integral.F.singularPts_y.length();

	DenseVector<Array<T> > AllSingularPoints_y(n_phi_y + n_F_y);

	std::merge(SingularPoints_phi_y.engine().data(),
			   SingularPoints_phi_y.engine().data() + n_phi_y,
		       integral.F.singularPts_y.engine().data(),
		       integral.F.singularPts_y.engine().data() + n_F_y,
		       AllSingularPoints_y.engine().data());

	T ret = 0.0;
	T left_x  = SingularPoints_phi_x(SingularPoints_phi_x.firstIndex());
	T right_x = SingularPoints_phi_x(SingularPoints_phi_x.lastIndex());
	T left_y  = SingularPoints_phi_y(SingularPoints_phi_y.firstIndex());
	T right_y = SingularPoints_phi_y(SingularPoints_phi_y.lastIndex());

	for (int i=AllSingularPoints_x.firstIndex(); i<AllSingularPoints_x.lastIndex(); ++i) {
		T a_x = AllSingularPoints_x(i), b_x = AllSingularPoints_x(i+1);
		if ((b_x <= left_x) || (a_x >= right_x)) continue;
		for (int j=AllSingularPoints_y.firstIndex(); j<AllSingularPoints_y.lastIndex(); ++j) {
			T a_y = AllSingularPoints_y(j), b_y = AllSingularPoints_y(j+1);
			if ((b_y <= left_y) || (a_y >= right_y)) continue;
			T tmp = integral.quadrature(a_x,b_x,a_y,b_y);
			ret += tmp;
			//std::cout << "[" << AllSingularPoints_x(i) << ", " << AllSingularPoints_x(i+1) << "] x "
			//          << "[" << AllSingularPoints_y(j) << ", " << AllSingularPoints_y(j+1) << "] : " << tmp << std::endl;
		}
	}

	return ret;


	return 0;
}

template <typename T, QuadratureType Quad, typename BasisFunc_x, typename BasisFunc_y>
Integral2D<T,Quad,BasisFunc_x,BasisFunc_y>::Integral2D(const Function2D<T> &_F,
										const BasisFunc_x &_phi_x, const BasisFunc_y &_phi_y)
    : F(_F), phi_x(_phi_x), phi_y(_phi_y), j_x(0), k_x(0), j_y(0), k_y(0), quadrature(*this)
{
}

template <typename T, QuadratureType Quad, typename BasisFunc_x, typename BasisFunc_y>
T
Integral2D<T,Quad,BasisFunc_x,BasisFunc_y>::operator()(int _j_x, int _k_x, int _j_y, int _k_y) const
{
	j_x = _j_x, k_x = _k_x, j_y = _j_y, k_y = _k_y;
	return _integrate(*this);
}

template <typename T, QuadratureType Quad, typename BasisFunc_x, typename BasisFunc_y>
T
Integral2D<T,Quad,BasisFunc_x,BasisFunc_y>::integrand(T x, T y) const
{
	return F(x,y) * phi_x(x,j_x,k_x) * phi_y(y,j_y,k_y);
}

}   //namespace lawa
