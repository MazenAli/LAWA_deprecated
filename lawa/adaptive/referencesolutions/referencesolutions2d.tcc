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
 *        HELMHOLTZ EXAMPLES 2D - Tensor examples
 * ================================================================================*/

template <typename T, typename Basis2D>
flens::DenseVector<Array<T> >
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::sing_pts_x;

template <typename T, typename Basis2D>
flens::DenseVector<Array<T> >
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::sing_pts_y;

template <typename T, typename Basis2D>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::deltas_x;

template <typename T, typename Basis2D>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::deltas_y;

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
    	if (nr==2) {
    		sing_pts_x.engine().resize(1); sing_pts_x(1) = 1./3.;
    		deltas_x.engine().resize(1,2); deltas_x(1,1) = 1./3.; deltas_x(1,2) = 0.2;
    	    sing_pts_y.engine().resize(1); sing_pts_y(1) = 1./3.;
    	    deltas_y.engine().resize(1,2); deltas_y(1,1) = 1./3.; deltas_y(1,2) = 1.;
    	}
    	if (nr==3) {
    		sing_pts_x.engine().resize(1); sing_pts_x(1) = 1./3.;
    		deltas_x.engine().resize(1,2); deltas_x(1,1) = 1./3.; deltas_x(1,2) = 4.;
    	}
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
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::exact_dx(T x, T y)
{
    return exact_x(x,1) * exact_y(y,0);
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::exact_dy(T x, T y)
{
	return exact_x(x,0) * exact_y(y,1);
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
        else if (nr==2) {
        	if (deriv_x==0)         return  std::exp(-0.1*fabs(x-1./3.));
        	else if (deriv_x==1) {
        	    if (x < 1./3.) 		return   0.1*std::exp(-0.1*fabs(x-1./3.));
        	    else		   		return  -0.1*std::exp(-0.1*fabs(x-1./3.));
        	}
        	else 			    	return   0.01*std::exp(-0.1*fabs(x-1./3.));
        }
        else if (nr==3) {
            if (deriv_x==0)         return  std::exp(-2.*fabs(x-1./3.));
            else if (deriv_x==1) {
            	if (x < 1./3.) 		return   2.*std::exp(-2.*fabs(x-1./3.));
            	else		   		return  -2.*std::exp(-2.*fabs(x-1./3.));
            }
            else {
            	if (x < 1./3.) 		return   4.*std::exp(-2.*fabs(x-1./3.));
            	else		   		return   4.*std::exp(-2.*fabs(x-1./3.));
            }
        }
    }
    else if ((domain1 == Periodic) && (domain2 == Interval)) {
        if (nr==1) {
            if (deriv_x==0)         return              std::cos(2*M_PI*x);
            else if (deriv_x==1)    return         -2*M_PI*std::sin(2*M_PI*x);
            else                    return -4*M_PI*M_PI*std::cos(2*M_PI*x);
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
            else                    return  std::exp(-0.5*(y-0.1)*(y-0.1)) * (4*0.5*0.5*(y-0.1)*(y-0.1)-2*0.5);
        }
        else if (nr==2) {
        	if (deriv_y==0)         return  std::exp(-0.5*fabs(y-1./3.));
        	else if (deriv_y==1) {
        	    if (y < 1./3.) return   0.5*std::exp(-0.5*fabs(y-1./3.));
        	    else		   return  -0.5*std::exp(-0.5*fabs(y-1./3.));
        	}
        	else			   return   0.25*std::exp(-0.5*fabs(y-1./3.));
        }
        else if (nr==3) {
        	if (deriv_y==0) 		return std::exp(-0.1*(y-1./3.)*(y-1./3.));
        	else if (deriv_y==1)	return -0.2*(y-1./3.)*std::exp(-0.1*(y-1./3.)*(y-1./3.));
        	else					return  (0.04*(y-1./3)*(y-1./3.)-0.2)*std::exp(-0.1*(y-1./3.)*(y-1./3.));
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
        if (nr==1) {
            //T a1 = -0.1;
            //T a2 = 0.1;
            T b1 = 0.1;
            T b2 = 0.5;
            T c1=1., c2= 1.;
            ret += c1*c1*std::sqrt(0.5*M_PI/(b1)) * c2*c2*std::sqrt(0.5*M_PI/(b2));
            ret += c1*c1*b1*b1*std::sqrt(0.5*M_PI)/std::pow(b1,1.5) * c2*c2*std::sqrt(0.5*M_PI/(b2));
            ret += c1*c1*std::sqrt(0.5*M_PI/(b1)) * c2*c2*b2*b2*std::sqrt(0.5*M_PI)/std::pow(b2,1.5);
            ret = sqrt(ret);
        }
        if (nr==2) {
            T i11   = 10.;
            T i22   = 2.;
            T ddi11 = 0.1;
            T ddi22 = 0.5;
            ret = sqrt(ddi11*i22 + i11*ddi22 + i11*i22);
        }
        if (nr==3) {
        	ret += 3.96332729760601*0.5 + 0.3963327297606012*0.5 + 3.96332729760601*2.;
        	return sqrt(ret);
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

/* ================================================================================
 *        HELMHOLTZ EXAMPLES 2D - Non tensor examples
 * ================================================================================*/

template <typename T, typename Basis2D>
flens::DenseVector<Array<T> >
ReferenceSolution2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::sing_pts_x;

template <typename T, typename Basis2D>
flens::DenseVector<Array<T> >
ReferenceSolution2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::sing_pts_y;

template <typename T, typename Basis2D>
T
ReferenceSolution2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::c;

template <typename T, typename Basis2D>
int
ReferenceSolution2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::nr;

template <typename T, typename Basis2D>
DomainType
ReferenceSolution2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::domain1;

template <typename T, typename Basis2D>
DomainType
ReferenceSolution2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::domain2;

template <typename T, typename Basis2D>
void
ReferenceSolution2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::setExample(int _nr,
                                                            const HelmholtzOperator2D<T,Basis2D> &a,
                                                            DomainType _domain1, DomainType _domain2)
{
    c=a.getc();
    assert(c>=0);
    nr=_nr;
    domain1 = _domain1;
    domain2 = _domain2;

    if ((domain1 == R) && (domain2 == R)) {
    	if (nr==2) {

    		sing_pts_x.engine().resize(9);
    		sing_pts_y.engine().resize(9);
    		sing_pts_x = -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5;
    		sing_pts_y = -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5;

/*
    		sing_pts_x.engine().resize(17);
    		sing_pts_y.engine().resize(17);
    		sing_pts_x = -0.3,-0.2, -0.1, 0., 0.099, 0.0999,  0.09999, 0.09999999, 0.1, 0.10000001, 0.10001, 0.1001, 0.101, 0.2, 0.3, 0.4, 0.5;
    		sing_pts_y = -0.3,-0.2, -0.1, 0., 0.099, 0.0999,  0.09999, 0.09999999, 0.1, 0.10000001, 0.10001, 0.1001, 0.101, 0.2, 0.3, 0.4, 0.5;
*/
    	}
    }
}

template <typename T, typename Basis2D>
T
ReferenceSolution2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::exact(T x, T y)
{
    return exact(x,y,0,0);
}

template <typename T, typename Basis2D>
T
ReferenceSolution2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::exact_dx(T x, T y)
{
    return exact(x,y,1,0);
}

template <typename T, typename Basis2D>
T
ReferenceSolution2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::exact_dy(T x, T y)
{
    return exact(x,y,0,1);
}

template <typename T, typename Basis2D>
T
ReferenceSolution2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::rhs(T x, T y)
{
    return 1.;
	//return -exact(x,y,2,0)-exact(x,y,0,2)+c*exact(x,y,0,0);
}

template <typename T, typename Basis2D>
T
ReferenceSolution2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::exact(T x, T y, int deriv_x, int deriv_y)
{
	if ((domain1=R) && (domain2==R)) {
		if (nr==1) {
			if ((deriv_x==0) && (deriv_y==0)) {
				return exp(-( 2*(x-0.1)*(x-0.1) + (x-0.1)*(y-0.1) + (y-0.1)*(y-0.1)  ) );
			}
			else if ((deriv_x==0) && (deriv_y==2)) {
				return (pow(-(x-0.1)-2*(y-0.1),2) - 2)*
				       exp(-( 2*(x-0.1)*(x-0.1) + (x-0.1)*(y-0.1) + (y-0.1)*(y-0.1)  ) );
			}
			else if ((deriv_x==2) && (deriv_y==0)) {
				return (pow(-4*(x-0.1)-(y-0.1),2) - 4)*
					   exp(-( 2*(x-0.1)*(x-0.1) + (x-0.1)*(y-0.1) + (y-0.1)*(y-0.1)  ) );
			}
			else {
				return 0.;
			}
		}
		else if (nr==2) {
			T XmA_p2_p_YmB_p2 = (x-0.1)*(x-0.1)+(y-0.1)*(y-0.1);
			if ((deriv_x==0) && (deriv_y==0)) {
				return exp(-sqrt( XmA_p2_p_YmB_p2 ));
			}
			else if ((deriv_x==0) && (deriv_y==1)) {
				return (-(y-0.1)*sqrt(XmA_p2_p_YmB_p2)/XmA_p2_p_YmB_p2)*exp(-sqrt( XmA_p2_p_YmB_p2 ));
			}
			else if ((deriv_x==1) && (deriv_y==0)) {
				return (-(x-0.1)*sqrt(XmA_p2_p_YmB_p2)/XmA_p2_p_YmB_p2)*exp(-sqrt( XmA_p2_p_YmB_p2 ));
			}
			else if ((deriv_x==0) && (deriv_y==2)) {
				return ( (y-0.1)*(y-0.1)/XmA_p2_p_YmB_p2 -
					     1./sqrt(XmA_p2_p_YmB_p2) + (y-0.1)*(y-0.1)/pow(XmA_p2_p_YmB_p2,1.5) )*
				       exp(-sqrt( XmA_p2_p_YmB_p2 ));
			}
			else if ((deriv_x==2) && (deriv_y==0)) {
				return ( (x-0.1)*(x-0.1)/XmA_p2_p_YmB_p2 -
					     1./sqrt(XmA_p2_p_YmB_p2) + (x-0.1)*(x-0.1)/pow(XmA_p2_p_YmB_p2,1.5) )*
					   exp(-sqrt( XmA_p2_p_YmB_p2 ));
			}
			else {
				return 0.;
			}
		}
		else {
			assert(0);
			return 0;
		}
	}
	else {
		assert(0);
		return 0;
	}
}

template <typename T, typename Basis2D>
T
ReferenceSolution2D<T,Basis2D,HelmholtzOperator2D<T,Basis2D> >::H1norm()
{
    T ret = 0.;
    if ((domain1==R) && (domain2==R)) {
        if (nr==1)             {
        	ret = 1.187410411723726 + 2.374820823447452 + 1.187410411723726;
        	ret = sqrt(ret);
        }
        else if (nr==2) {
        	ret = 1.570795505591071 + 0.7853981678060725 + 0.7853981678060725;
        	ret = sqrt(ret);
        }
    }
    return ret;
}

/* ================================================================================
 *        SPACE-TIME HEAT EXAMPLES 2D
 * ================================================================================*/
 
template <typename T, typename Basis2D>
flens::DenseVector<Array<T> >
ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::sing_pts_t;

template <typename T, typename Basis2D>
flens::DenseVector<Array<T> >
ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::sing_pts_x;

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::c;

template <typename T, typename Basis2D>
int
ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::nr;

template <typename T, typename Basis2D>
DomainType
ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::domain1;

template <typename T, typename Basis2D>
DomainType
ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::domain2;

template <typename T, typename Basis2D>
void
ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::setExample(int _nr,
															const SpaceTimeHeatOperator1D<T,Basis2D> &a,
															DomainType _domain1, DomainType _domain2)
{
    c=a.getc();
    assert(c>=0);
    nr=_nr;
    domain1 = _domain1;
    domain2 = _domain2;

    if ((domain1 == Periodic) && (domain2 == Interval)) {
        switch(nr){
            case 1:
                sing_pts_x.engine().resize(2);
                sing_pts_x = 0., 1.;
                break;
            case 2: 
                sing_pts_t.engine().resize(3);
                sing_pts_t = 0., 0.5, 1.;
            case 3: 
                sing_pts_t.engine().resize(3);
                sing_pts_t = 0., 2.-std::sqrt(2), 1.;
            default: break;
        }
    }
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::exact(T t, T x)
{
	return exact_t(t,0)*exact_x(x,0);
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::dx_exact(T t, T x)
{
    return exact_t(t,0) * exact_x(x,1);
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::_exact_t(T t)
{
	return exact_t(t,0);
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::exact_t(T t)
{
	return exact_t(t,0);
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::exact_x(T x)
{
	return exact_x(x,0);
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::rhs_t(T t)
{
	return exact_t(t,1);
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::rhs_x(T x)
{
	return -c*exact_x(x,2);
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::dd_exact_x(T x)
{
    return exact_x(x,2);
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::exact_t(T t, int deriv_t)
{
    if ((domain1 == Periodic) && (domain2 == Interval)) {
		switch(nr){
		    case 1:
    			if (deriv_t==0) 		return  			std::cos(2*M_PI*t);
    			else if (deriv_t==1)	return 		-2*M_PI*std::sin(2*M_PI*t);
                break;
            case 2:
                if (deriv_t==0){
                    if(t < 0.5) return  t + 0.25;
                    else        return -t + 1.25;
                }
    			else if (deriv_t==1){
    			    if(t < 0.5) return  1;
                    else        return -1;
    			}
                break;
            case 3:
                if (deriv_t==0){
                    if(t < 2. - std::sqrt(2.)) return  t*t + 0.5;
                    else        return  2*(t-1)*(t-1) + 0.5;
                }
    			else if (deriv_t==1){
    			    if(t < 2.-std::sqrt(2)) return  2*t;
                    else        return 4*(t-1);
    			}
                break;
            default: std::cerr << "Example does not exist!" << std::endl; exit(1);
		}
	}
	else { std::cerr << "Example does not exist!" << std::endl; exit(1); }
    return 0;
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::exact_x(T x, int deriv_x)
{
    if ((domain1 == Periodic) && (domain2 == Interval)) {
		switch(nr){
		    case 1:
    			if (deriv_x==0) 		return  -4*(x-0.5)*(x-0.5)+1;
    			else if (deriv_x==1)	return	-8*(x-0.5);
    			else 					return  -8;
                break;
            case 2:
                if (deriv_x==0) 		return  8*std::pow(x-0.5, 3) - 2*x*x +1;
    			else if (deriv_x==1)	return	24*(x-0.5)*(x-0.5) -4*x;
    			else 					return  48*(x-0.5) - 4;
                break;
            case 3:
                if (deriv_x==0) 		return  8*std::pow(x-0.5, 3) - 2*x*x +1;
    			else if (deriv_x==1)	return	24*(x-0.5)*(x-0.5) -4*x;
    			else 					return  48*(x-0.5) - 4;
                break;
            default: std::cerr << "Example does not exist!" << std::endl; exit(1); 
		}
	}
	else { std::cerr << "Example does not exist!" << std::endl; exit(1); }
    return 0;
}

template <typename T, typename Basis2D>
T
ReferenceSolutionTensor2D<T,Basis2D,SpaceTimeHeatOperator1D<T,Basis2D> >::H1_t_norm(T t){
    T ret = 0.0;
    if ((domain1 == Periodic) && (domain2 == Interval)) {
		switch(nr){
		    case 1:
                ret = exact_t(t) * sqrt(5.866666666);
                break;
            case 2:
                break;
            default: break; 
		}
	}
	return ret;
}

}	//namespace lawa
