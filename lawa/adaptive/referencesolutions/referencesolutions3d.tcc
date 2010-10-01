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

/* ================================================================================
 *        HELMHOLTZ EXAMPLES 3D - Tensor examples
 * ================================================================================*/

template <typename T, typename Basis3D>
flens::DenseVector<Array<T> >
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::sing_pts_x;

template <typename T, typename Basis3D>
flens::DenseVector<Array<T> >
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::sing_pts_y;

template <typename T, typename Basis3D>
flens::DenseVector<Array<T> >
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::sing_pts_z;

template <typename T, typename Basis3D>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::deltas_x;

template <typename T, typename Basis3D>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::deltas_y;

template <typename T, typename Basis3D>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::deltas_z;

template <typename T, typename Basis3D>
T
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::c;

template <typename T, typename Basis3D>
int
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::nr;

template <typename T, typename Basis3D>
DomainType
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::domain1;

template <typename T, typename Basis3D>
DomainType
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::domain2;

template <typename T, typename Basis3D>
DomainType
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::domain3;

template <typename T, typename Basis3D>
void
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::setExample(int _nr,
                                                            const HelmholtzOperator3D<T,Basis3D> &a,
                                                            DomainType _domain1, DomainType _domain2,
                                                            DomainType _domain3)
{
    c=a.getc();
    assert(c>=0);
    nr=_nr;
    domain1 = _domain1;
    domain2 = _domain2;
    domain3 = _domain3;


}

template <typename T, typename Basis3D>
T
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::exact(T x, T y, T z)
{
    return exact_x(x,0)*exact_y(y,0)*exact_z(z,0);
}

template <typename T, typename Basis3D>
T
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::exact_dx(T x, T y, T z)
{
    return exact_x(x,1) * exact_y(y,0) * exact_z(z,0);
}

template <typename T, typename Basis3D>
T
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::exact_dy(T x, T y, T z)
{
	return exact_x(x,0) * exact_y(y,1) * exact_z(z,0);
}

template <typename T, typename Basis3D>
T
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::exact_dz(T x, T y, T z)
{
	return exact_x(x,0) * exact_y(y,0) * exact_z(z,1);
}

template <typename T, typename Basis3D>
T
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::exact_x(T x)
{
    return exact_x(x,0);
}

template <typename T, typename Basis3D>
T
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::exact_y(T y)
{
    return exact_y(y,0);
}

template <typename T, typename Basis3D>
T
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::exact_z(T z)
{
    return exact_z(z,0);
}

template <typename T, typename Basis3D>
T
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::rhs_x(T x)
{
    return -exact_x(x,2) + (1./3.)*c*exact_x(x,0);
}

template <typename T, typename Basis3D>
T
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::rhs_y(T y)
{
    return -exact_y(y,2) + (1./3.)*c*exact_y(y,0);
}

template <typename T, typename Basis3D>
T
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::rhs_z(T z)
{
    return -exact_z(z,2) + (1./3.)*c*exact_z(z,0);
}

template <typename T, typename Basis3D>
T
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::exact_x(T x, int deriv_x)
{
    if ((domain1==R) && (domain2==R) && (domain3==R)) {
        if (nr==1) {
            if (deriv_x==0)         return  std::exp(-0.1*(x+0.1)*(x+0.1));
            else if (deriv_x==1)    return -std::exp(-0.1*(x+0.1)*(x+0.1)) * 2*0.1*(x+0.1);
            else                     return std::exp(-0.1*(x+0.1)*(x+0.1)) * (4*0.1*0.1*(x+0.1)*(x+0.1)-2*0.1);
        }
    }
    return 0;
}

template <typename T, typename Basis3D>
T
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::exact_y(T y, int deriv_y)
{
    if ((domain1==R) && (domain2==R) && (domain3==R)) {
        if (nr==1) {
            if (deriv_y==0)         return  std::exp(-0.5*(y-0.1)*(y-0.1));
            else if (deriv_y==1)    return -std::exp(-0.5*(y-0.1)*(y-0.1)) * 2*0.5*(y-0.1);
            else                    return  std::exp(-0.5*(y-0.1)*(y-0.1)) * (4*0.5*0.5*(y-0.1)*(y-0.1)-2*0.5);
        }
    }
    return 0;
}

template <typename T, typename Basis3D>
T
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::exact_z(T z, int deriv_z)
{
    if ((domain1==R) && (domain2==R) && (domain3==R)) {
        if (nr==1) {
            if (deriv_z==0)         return  std::exp(-(z-0.1)*(z-0.1));
            else if (deriv_z==1)    return -std::exp(-(z-0.1)*(z-0.1)) * 2*(z-0.1);
            else                    return  std::exp(-(z-0.1)*(z-0.1)) * (4*(z-0.1)*(z-0.1)-2);
        }
    }
    return 0;
}

template <typename T, typename Basis3D>
T
ReferenceSolutionTensor3D<T,Basis3D,HelmholtzOperator3D<T,Basis3D> >::H1norm()
{
    T ret = 0.;
    if ((domain1==R) && (domain2==R) && (domain3==R)) {
        if (nr==1) {
            T i11   = 3.963327297606011;
            T i22   = 1.772453850905516;
            T i33   = 1.2533141373155;
            T ddi11 = 0.3963327297606013;
            T ddi22 = 0.886226925452758;
            T ddi33 = 1.2533141373155;

            ret = sqrt(ddi11*i22*i33 + i11*ddi22*i33 + i11*i22*ddi33 + i11*i22*i33);
        }
    }
    return ret;
}

}	//namespace lawa
