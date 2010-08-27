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

template <typename T, typename Basis2D>
flens::DenseVector<Array<T> >
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::sing_pts_x;

template <typename T, typename Basis2D>
flens::DenseVector<Array<T> >
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::sing_pts_y;

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::c;

template <typename T, typename Basis2D>
int
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::nr;

template <typename T, typename Basis2D>
DomainType
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::domain1;

template <typename T, typename Basis2D>
DomainType
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::domain2;

template <typename T, typename Basis2D>
void
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::setExample(int _nr,
                                                            const HelmholtzOperator2D<T,Basis2D> &a,
                                                            DomainType _domain1, DomainType _domain2)
{
    c=a.getc();
    assert(c>=0);
    nr=_nr;
    domain1 = _domain1;
    domain2 = _domain2;

    if ((domain1 == R) && (domain2 == R)) {

    }
    else if ((domain1 == Periodic) && (domain2 == Interval)) {

    }
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::exact(T x, T y)
{
    return exact_x(x,0)*exact_y(y,0);
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::exact_x(T x)
{
    return exact_x(x,0);
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::exact_y(T y)
{
    return exact_y(y,0);
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::rhs_x(T x)
{
    return -exact_x(x,2) + 0.5*c*exact_x(x,0);
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::rhs_y(T y)
{
    return -exact_y(y,2) + 0.5*c*exact_y(y,0);
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::exact_x(T x, int deriv_x)
{
    if ((domain1==R) && (domain2==R)) {
        if (nr==1) {
            if (deriv_x==0)         return  std::exp(-0.1*(x+0.1)*(x+0.1));
            else if (deriv_x==1)    return -std::exp(-0.1*(x+0.1)*(x+0.1)) * 2*0.1*(x+0.1);
            else                     return  std::exp(-0.1*(x+0.1)*(x+0.1)) * (4*0.1*0.1*(x+0.1)*(x+0.1)-2*0.1);
        }
    }
    else if ((domain1 == Periodic) && (domain2 == Interval)) {
        if (nr==1) {
            if (deriv_x==0)         return              std::cos(2*M_PI*x);
            else if (deriv_x==1)    return         -2*M_PI*std::sin(2*M_PI*x);
            else                     return -4*M_PI*M_PI*std::cos(2*M_PI*x);
        }
        else { std::cerr << "Example does not exist!" << std::endl; exit(1); }
        return 0;
    }
    else { std::cerr << "Example does not exist!" << std::endl; exit(1); }
    return 0;
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::exact_y(T y, int deriv_y)
{
    if ((domain1==R) && (domain2==R)) {
        if (nr==1) {
            if (deriv_y==0)         return  std::exp(-0.5*(y-0.1)*(y-0.1));
            else if (deriv_y==1)    return -std::exp(-0.5*(y-0.1)*(y-0.1)) * 2*0.5*(y-0.1);
            else                     return  std::exp(-0.5*(y-0.1)*(y-0.1)) * (4*0.5*0.5*(y-0.1)*(y-0.1)-2*0.5);
        }
    }
    else if ((domain1 == Periodic) && (domain2 == Interval)) {
        if (nr==1) {
            if (deriv_y==0)         return  -4*(y-0.5)*(y-0.5)+1;
            else if (deriv_y==1)    return    -8*(y-0.5);
            else                     return  -8;
        }
        else { std::cerr << "Example does not exist!" << std::endl; exit(1); }
        return 0; 
    }
    else { std::cerr << "Example does not exist!" << std::endl; exit(1); }
    return 0;
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::H1norm()
{
    T ret = 0.;
    if ((domain1==R) && (domain2==R)) {
        if (nr==1)             {
            T a1 = -0.1;
            T a2 = 0.1;
            T b1 = 0.1;
            T b2 = 0.5;
            T c1=1., c2= 1.;
            ret += c1*c1*std::sqrt(0.5*M_PI/(b1)) * c2*c2*std::sqrt(0.5*M_PI/(b2));
            ret += c1*c1*b1*b1*std::sqrt(0.5*M_PI)/std::pow(b1,1.5) * c2*c2*std::sqrt(0.5*M_PI/(b2));
            ret += c1*c1*std::sqrt(0.5*M_PI/(b1)) * c2*c2*b2*b2*std::sqrt(0.5*M_PI)/std::pow(b2,1.5);
            ret = sqrt(ret);
        }
    }
    else if ((domain1 == Periodic) && (domain2 == Interval)) {
        if (nr==1) {
            ret +=  0.2666666666666666;
            ret += 10.527578027828639;
            ret +=  2.666666666666666;
            ret = sqrt(ret);
        }
    }
    return ret;
}

}    //namespace lawa
