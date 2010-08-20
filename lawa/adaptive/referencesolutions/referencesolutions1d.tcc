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

template <typename T, typename Basis>
flens::DenseVector<Array<T> >
ReferenceSolution1D<T,Basis,HelmholtzOperator1D<T,Basis> >::sing_pts;

template <typename T, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
ReferenceSolution1D<T,Basis,HelmholtzOperator1D<T,Basis> >::deltas;

template <typename T, typename Basis>
T
ReferenceSolution1D<T,Basis,HelmholtzOperator1D<T,Basis> >::c;

template <typename T, typename Basis>
int
ReferenceSolution1D<T,Basis,HelmholtzOperator1D<T,Basis> >::nr;

template <typename T, typename Basis>
DomainType
ReferenceSolution1D<T,Basis,HelmholtzOperator1D<T,Basis> >::domain;


template <typename T, typename Basis>
void
ReferenceSolution1D<T,Basis,HelmholtzOperator1D<T,Basis> >::setExample(int _nr,
						const HelmholtzOperator1D<T,Basis> &a, DomainType _domain)
{
    c=a.getc();
    assert(c>=0);
    nr=_nr;
    domain = _domain;

    if (domain == Interval) {

    }
    else if (domain == Periodic) {

    }
    else { 	//domain == R
    	if (nr==1) {
    	//	deltas.engine().resize(1,2);
    	//	deltas(1,1) = 0.01; deltas(1,2) = 0.;
    	}
    	else if (nr==2) {
    		sing_pts.engine().resize(1);
    		sing_pts(1) = 0.01;
    		deltas.engine().resize(1,2);
    		deltas(1,1) = 0.01; deltas(1,2) = 2.;
    	}
    }
}

template <typename T, typename Basis>
T
ReferenceSolution1D<T,Basis,HelmholtzOperator1D<T,Basis> >::exact(T x, int deriv)
{
	if (domain==Interval) {
		if (nr==1) {
			T a = 100;
			if (deriv==0) 		return exp(-a*(x-0.5)*(x-0.5));
			else if (deriv==1)	return -exp(-a*(x-0.5)*(x-0.5))*(2*a*(x-0.5));
			else if (deriv==2)	return exp(-a*(x-0.5)*(x-0.5)) * ((2*a*(x-0.5))*(2*a*(x-0.5))-2*a);
			else				assert(0);
		}
		else 					assert(0);
	}
	else if (domain==Periodic) {
		if (nr==1) {
			if (deriv == 0)		 {	return std::cos(2*M_PI*x);		}
			else if (deriv == 1) {	return -2*M_PI*std::sin(2*M_PI*x); }
			else if (deriv == 2) {	return -4*M_PI*M_PI*std::cos(2*M_PI*x);  }
		}
		else { std::cout << "ReferenceSolution<1d> exits!" << std::endl; exit(1); }
	}
	else { //domain == R
		if (nr==1) {
			if (deriv == 0)		 {	return 10.*std::exp(-0.1*(x-0.1)*(x-0.1));		}
			else if (deriv == 1) {	return -2*(x-0.1)*std::exp(-0.1*(x-0.1)*(x-0.1)); }
			else if (deriv == 2) {	return (4*0.1*(x-0.1)*(x-0.1)-2)*std::exp(-0.1*(x-0.1)*(x-0.1));  }
		}
		else if (nr==2) {
			if (deriv == 0)		 {	return 10.*std::exp(-0.1*fabs(x-0.01));		}
			else if (deriv == 1) {
				if (x > 0.01) 		return -std::exp(-0.1*(x-0.01));
				else	  			return  std::exp( 0.1*(x-0.01));
			}
			else if (deriv == 2) {	//only piecewise, delta is required!!!
				if (x > 0.01) 		return 0.1*std::exp(-0.1*(x-0.01));
				else	  			return 0.1*std::exp( 0.1*(x-0.01));
			}
		}
		else { std::cout << "ReferenceSolution<1d> exits!" << std::endl; exit(1); }
	}
}

template <typename T, typename Basis>
T
ReferenceSolution1D<T,Basis,HelmholtzOperator1D<T,Basis> >::exact(T x)
{
    return ReferenceSolution1D<T,Basis,HelmholtzOperator1D<T,Basis> >::exact(x, 0);
}

template <typename T, typename Basis>
T
ReferenceSolution1D<T,Basis,HelmholtzOperator1D<T,Basis> >::rhs(T x)
{
    return -exact(x,2) + c * exact(x,0);
}

template <typename T, typename Basis>
T
ReferenceSolution1D<T,Basis,HelmholtzOperator1D<T,Basis> >::H1norm()
{
	if (domain==R) {
		if (nr==1) 			{ return std::sqrt(100*std::sqrt(0.5*M_PI/(0.1)) + std::sqrt(0.5*M_PI)/std::pow(0.1,(T)1.5) );	}
		else if (nr ==2 )	{ return std::sqrt(1000. + 10.); }
		else {	std::cerr << "ReferenceSolution<1d> exits!" << std::endl; exit(1);	}
	}
	else {
		return 0.;
	}
}

} //namespace lawa
