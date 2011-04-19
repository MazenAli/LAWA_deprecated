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

template <typename T>
int
PDEReferenceSolutions1D<T>::nr;

template <typename T>
T
PDEReferenceSolutions1D<T>::diffusion;

template <typename T>
T
PDEReferenceSolutions1D<T>::convection;

template <typename T>
T
PDEReferenceSolutions1D<T>::reaction;


template <typename T>
DomainType
PDEReferenceSolutions1D<T>::domain;

template <typename T>
flens::DenseVector<Array<T> >
PDEReferenceSolutions1D<T>::sing_pts;

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
PDEReferenceSolutions1D<T>::deltas;


template <typename T>
void
PDEReferenceSolutions1D<T>::setExample(int _nr, T _diffusion, T _convection, T _reaction,
                                       DomainType _domain)
{
    nr=_nr;
    diffusion = _diffusion;
    convection = _convection;
    reaction = _reaction;
    domain = _domain;

    if (domain == Interval) {

    }
    else if (domain == Periodic) {

    }
    else {     //domain == R
        assert(nr>=1);
        assert(nr<=8);
        if (nr==1) {

        }
        else if (nr==2) {
            sing_pts.engine().resize(1);
            sing_pts(1) = 0.01;
            deltas.engine().resize(1,2);
            deltas(1,1) = 0.01; deltas(1,2) = diffusion*2.;
        }
        else if (nr==3) {
            sing_pts.engine().resize(2);
            sing_pts(1) = -M_PI/8.;
            sing_pts(2) =  M_PI/8.;
            deltas.engine().resize(2,2);
            deltas(1,1) = -M_PI/8.; deltas(1,2) = -diffusion*4.;
            deltas(2,1) =  M_PI/8.; deltas(2,2) = diffusion*4.;
        }
        else if (nr==4) {    //second derivative not continuous on R!!
            sing_pts.engine().resize(2);
            sing_pts(1) = -0.4;
            sing_pts(2) =  0.9;
        }
        else if (nr==5) {    //second derivative not continuous on R!!
            sing_pts.engine().resize(2);
            sing_pts(1) = -0.4;
            sing_pts(2) =  0.9;
            deltas.engine().resize(2,2);
			deltas(1,1) = -0.4; deltas(1,2) = diffusion*20.8;
			deltas(2,1) =  0.9; deltas(2,2) = diffusion*105.3;
        }
        else if (nr==6) {
        	sing_pts.engine().resize(1);
        	sing_pts(1) = 0.0001;
        	deltas.engine().resize(1,2);
        	//deltas(1,1) = 0.0001; deltas(1,2) = 3.;
        	deltas(1,1) = 0.0001; deltas(1,2) = diffusion*9.;
        }
        else if (nr==7) {
            sing_pts.engine().resize(1);
            sing_pts(1) = 0.0001;
            deltas.engine().resize(1,2);
            deltas(1,1) = 0.0001; deltas(1,2) = diffusion*3.;
        }
        else if (nr==8) {
            sing_pts.engine().resize(2);
            sing_pts(1) = -16./3.;
            sing_pts(2) =  16./3.;
            deltas.engine().resize(2,2);
            deltas(1,1) = -16./3.; deltas(1,2) = 22773./6391.;
            deltas(2,1) =  16./3.; deltas(2,2) = 16382./6391.;
        }
    }
}

template <typename T>
T
PDEReferenceSolutions1D<T>::exact(T x, int deriv)
{
    if (domain==Interval) {
        if (nr==1) {
            T a = 100;
            if (deriv==0)         return exp(-a*(x-0.5)*(x-0.5));
            else if (deriv==1)    return -exp(-a*(x-0.5)*(x-0.5))*(2*a*(x-0.5));
            else if (deriv==2)    return exp(-a*(x-0.5)*(x-0.5))*((2*a*(x-0.5))*(2*a*(x-0.5))-2*a);
            else                  assert(0);
        }
        else                     assert(0);
    }
    else if (domain==Periodic) {
        if (nr==1) {
            if (deriv == 0)         return std::cos(2*M_PI*x);
            else if (deriv == 1)    return -2*M_PI*std::sin(2*M_PI*x);
            else if (deriv == 2)    return -4*M_PI*M_PI*std::cos(2*M_PI*x);
            else                    assert(0);
        }
        else                         assert(0);
    }
    else { //domain == R
        if 	    (nr==1) {
            if (deriv == 0)         return 10.*std::exp(-0.1*(x-0.1)*(x-0.1));
            else if (deriv == 1)    return -2*(x-0.1)*std::exp(-0.1*(x-0.1)*(x-0.1));
            else if (deriv == 2)    return (4*0.1*(x-0.1)*(x-0.1)-2)*std::exp(-0.1*(x-0.1)*(x-0.1));
            else                    assert(0);
        }
        else if (nr==2) {
            if (deriv == 0)             return 10.*std::exp(-0.1*fabs(x-0.01));
            else if (deriv == 1) {
                if (x > 0.01)         return -std::exp(-0.1*(x-0.01));
                else                  return  std::exp( 0.1*(x-0.01));
            }
            else if (deriv == 2) {    //only piecewise, delta is required!!!
                if (x > 0.01)         return 0.1*std::exp(-0.1*(x-0.01));
                else                  return 0.1*std::exp( 0.1*(x-0.01));
            }
            else                    assert(0);
        }
        else if (nr == 3){
            if (deriv == 0) {
                if (x<-M_PI/8.)    	     return -1./std::pow(x+M_PI/8.-1,(T)4);
                else if (x>M_PI/8.)      return  1./std::pow(x-M_PI/8.+1,(T)4);
                else                	 return std::sin(4*x);
             }
             else if (deriv == 1) {
                 if (x<-M_PI/8.)    	 return  4./std::pow(x+M_PI/8.-1,(T)5);
                 else if (x>M_PI/8.)     return -4./std::pow(x-M_PI/8.+1,(T)5);
                 else                	 return 4*std::cos(4*x);
             }
             else if (deriv == 2) {
                 if (x<-M_PI/8.)    	 return -20./std::pow(x+M_PI/8.-1,(T)6);
                 else if (x>M_PI/8.)     return  20./std::pow(x-M_PI/8.+1,(T)6);
                 else                	 return -16*std::sin(4*x);

             }
             else                         assert(0);
        }
        else if (nr==4) {
            if (deriv==0) {
                if ((x>-0.4) && (x<0.9))  return 100*(x+0.4)*(x+0.4)*(x-0.9)*(x-0.9);
                else                       return 0;
            }
            else if (deriv==1) {
                if ((x>-0.4) && (x<0.9))  return   200*(x+0.4)*(x-0.9)*(x-0.9)
                                                 + 200*(x+0.4)*(x+0.4)*(x-0.9);
                else                       return 0;
            }
            else if (deriv==2) {
                if ((x>-0.4) && (x<0.9))  return   200*(x-0.9)*(x-0.9) + 800.*(x+0.4)*(x-0.9)
                                                 + 200.*(x+0.4)*(x+0.4);
                else                       return 0;
            }
            else                         assert(0);
        }
        else if (nr==5) {
            if (deriv==0) {
            	if ((x>-0.4) && (x<0.9))   return 100*(x+0.4)*x*x*(x-0.9);
            	else return 0;
            }
            else if (deriv==1) {
                if ((x>-0.4) && (x<0.9))   return   100*(x-0.9)*x*x + 200*(x-0.9)*x*(x+0.4)
                                                  + 100*x*x*(x+0.4);
                else return 0;
            }
            else if (deriv==2) {
                if ((x>-0.4) && (x<0.9))   return   400*x*(x-0.9) + 200*x*x + 200*(x-0.9)*(x+0.4)
                                                  + 400*(x+0.4)*x;
                else                       return 0;
            }
            else                         assert(0);
        }
        else if (nr==6) {
        	if (deriv==0) {
        		if (x<=0.0001)     	return 1./(std::pow(x-0.0001-1,(T)4));
                else 			 	return 1./(std::pow(x-0.0001+1,(T)5));
            }
            else if (deriv==1) {
            	if (x<=0.0001)      return -4./(std::pow(x-0.0001-1,(T)5));
                else 		 	 	return -5./(std::pow(x-0.0001+1,(T)6));
            }
            else if (deriv==2) {
            	if (x<=0.0001)      return 20./(std::pow(x-0.0001-1,(T)6));
            	else 		 	 	return 30./(std::pow(x-0.0001+1,(T)7));
            }
        }
        else if (nr==7) {
            T tmp = x-0.0001-1;
            if (deriv==0) {
                if (x<=0.0001)      return exp(2*(x-0.0001));
                else                return exp(-(x-0.0001));
            }
            else if (deriv==1) {
                if (x<=0.0001)      return 2*exp(2*(x-0.0001));
                else                return -exp(-(x-0.0001));
            }
            else if (deriv==2) {
                if (x<=0.0001)      return 4*exp(2*(x-0.0001));
                else                return exp(-(x-0.0001));
            }
        }
        else if (nr == 8){
            if (deriv == 0) {
                T coeff1 = 675./12782., coeff2 = -3209./6391.;
                if (x<-16./3.)         return -1./std::pow(x-1+16./3.,(T)3.);
                else if (x>16/3.)      return  1./std::pow(x+1-16./3.,(T)2.);
                else                 return coeff1*x*x+coeff2;
            }
            else if (deriv == 1) {
                T coeff1 = 675./12782.;
                if (x<-16./3.)             return  3./std::pow(x-1+16./3.,(T)4.);
                else if (x>16./3.)      return -2./std::pow(x+1-16./3.,(T)3.);
                else                     return 2.*coeff1*x;
            }
            else if (deriv == 2) {
                T coeff1 = 675./12782.;
                if (x<-16./3.)             return  -12./std::pow(x-1+16./3.,(T)5.);
                else if (x>16./3.)      return 6./std::pow(x+1-16./3.,(T)4.);
                else                     return 2.*coeff1;
            }
            else                         assert(0);
        }
        else {
            assert(0);
            return 0;
        }
    }
    assert(0);
    return 0;
}

template <typename T>
T
PDEReferenceSolutions1D<T>::exact(T x)
{
    return PDEReferenceSolutions1D<T>::exact(x, 0);
}

template <typename T>
T
PDEReferenceSolutions1D<T>::d_exact(T x)
{
    return PDEReferenceSolutions1D<T>::exact(x, 1);
}

template <typename T>
T
PDEReferenceSolutions1D<T>::rhs(T x)
{
	return -diffusion*exact(x,2) + convection*exact(x,1) + reaction*exact(x,0);
}

template <typename T>
T
PDEReferenceSolutions1D<T>::H1norm()
{
    if (domain==Interval) {
        if (nr==1)          return sqrt(0.12533141373155 + 12.533141373155);
        else                return 0.;
    }
    else if (domain==Periodic) {
        if (nr==1)            { return sqrt(0.5 + 2.*M_PI*M_PI); }
        else                return 0.;
    }
    else  {    // (domain==R)
        if (nr==1)          return   std::sqrt(100*std::sqrt(0.5*M_PI/(0.1))
                                   + std::sqrt(0.5*M_PI)/std::pow(0.1,(T)1.5) );
        else if (nr==2)     return std::sqrt(1000. + 10.);
        else if (nr==3)     return std::sqrt( (16+7.*M_PI)/56. + (2./9.)*(16.+9.*M_PI)   );
//      else if (nr==3)     return std::sqrt(1637492008./612673215.  + 656353759./204224405.);
        else if (nr==4)     return std::sqrt(168.3253868730168 + 1195.209847619049);
        else if (nr==5)		return std::sqrt(43.3676117539685 + 980.929114285711);
//      else if (nr==6)		return std::sqrt(0.75 + 1.5);
        else if (nr==6)		return std::sqrt(0.253968253968254 + 4.05050505050505);
//      else if (nr==6)		return std::sqrt(0.5333333333333333 + 2.085714285714286);
    }
    assert(0);
    return 0;
}

} //namespace lawa