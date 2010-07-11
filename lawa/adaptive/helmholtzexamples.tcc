/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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
DenseVector<Array<T> >
HelmholtzExamples<T>::u_sing_pts;

template <typename T>
T
HelmholtzExamples<T>::alpha;

template <typename T>
T
HelmholtzExamples<T>::beta;

template <typename T>
T
HelmholtzExamples<T>::energyNorm;

template <typename T>
T
HelmholtzExamples<T>::energyFunctional;

template <typename T>
int
HelmholtzExamples<T>::nr;

template <typename T>
void
HelmholtzExamples<T>::setExample(int _nr, T _alpha, T _beta)
{
    alpha=_alpha;
    beta=_beta;
    assert(beta>=0);
    assert(beta+0.5*alpha>=0);
    nr=_nr;
    assert(1<=nr && nr <=4);

    if (nr==1) {
        energyNorm=-100*exp(-50)+5*sqrt(2*M_PI)*erf(5*sqrt(2))
            + alpha * 0
            + beta * 1.0/20.0 * sqrt(2*M_PI)*erf(5*sqrt(2));
    } else if (nr==2) {
        energyNorm=40.0/3.0*(std::exp(5)+1)/(std::exp(5)-1)
            + alpha * 0
            + beta * (-4.0/15.0
                 * (-1 + 60*std::exp(10) + 8*std::exp(5) + std::exp(20) - 8*std::exp(15))
                  /(-std::exp(20) + 4*std::exp(15) - 6*std::exp(10) + 4*std::exp(5) - 1)
                     );
    } else if (nr==3) {
        energyNorm=2.88
            + alpha * 0
            + beta * 81.0/3125.0;
        u_sing_pts.engine().resize(2);
        u_sing_pts(1) = 0.3;
        u_sing_pts(2) = 0.7;
    } else if (nr==4) {
        energyNorm=-0.5*100.0*std::exp(100*(-0.6))
                        + 100.0
                        - 0.5*100.0*std::exp(100*(-1.4))
            + alpha * 0
            + beta * (-0.5/100.0*std::exp(100*(-0.6))
                        + 1.0/100.0
                        - 0.5/100.0*std::exp(100*(-1.4)));
        u_sing_pts.engine().resize(1);
        u_sing_pts(1) = 0.3;
    }
    energyFunctional=-0.5*energyNorm;
}

template <typename T>
T
HelmholtzExamples<T>::exact(T x)
{
    return HelmholtzExamples<T>::exact(x, 0);
}

template <typename T>
T
HelmholtzExamples<T>::exact(T x, int deriv)
{
    if (nr==1) {
        T a=100;
        if (deriv==0) {
            return exp(-a*(x-0.5)*(x-0.5));
        } else if (deriv==1) {
            return -exp(-a*(x-0.5)*(x-0.5))*(2*a*(x-0.5));
        } else if (deriv==2) {
            return exp(-a*(x-0.5)*(x-0.5)) * ((2*a*(x-0.5))*(2*a*(x-0.5))-2*a);
        } else {
            assert(0);
        }
    } else if (nr==2) {
        T a = 5.0;
        if (deriv==0) {
            return 4 * (exp(a*x)-1)/(exp(a)-1) * (1 - (exp(a*x)-1)/(exp(a)-1));
        } else if (deriv==1) {
            return (4*a*exp(a*x)*(exp(a)-2*exp(a*x)+1)) / ((exp(a)-1)*(exp(a)-1));
        } else if (deriv==2) {
            return (4 * a*a * exp(a*x) * (exp(a) - 4*exp(a*x) + 1) ) / ( (exp(a)-1)*(exp(a)-1) );
        } else {
            assert(0);
        }
    } else if (nr==3) {
        if (deriv==0) {
            if (x<=0.3) {
                return 4*x*x;
            } else if ((x>0.3) && (x<=0.7)) {
                return 9*(x-0.5)*(x-0.5);
            } else {
                return 4*(1-x)*(1-x);
            }
        } else if (deriv==1) {
            if (x<=0.3) {
                return 8*x;
            } else if ((x>0.3) && (x<=0.7)) {
                return 18*(x-0.5);
            } else {
                return -8*(1-x);
            }
        } else if (deriv==2) {
            if (x<=0.3) {
                return 8;
            } else if ((x>0.3) && (x<=0.7)) {
                return 18;
            } else {
                return 8;
            }
        } else {
            assert(0);
        }
    } else if (nr==4) {
        if (deriv==0) {
            if (x<=0.3) {
                return exp(100*(x-0.3));
            } else {
                return exp(100*(0.3-x));
            }
        } else if (deriv==1) {
            if (x<=0.3) {
                return  100*exp(100*(x-0.3));
            } else {
                return -100*exp(100*(0.3-x));
            }
        } else if (deriv==2) {
            if (x<=0.3) {
                return 100*100*exp(100*(x-0.3));
            } else {
                return 100*100*exp(100*(0.3-x));
            }
        } else {
            assert(0);
        }
    }

    assert(0);
    return 0;
}

template <typename T>
T
HelmholtzExamples<T>::rhs(T x)
{
    return -exact(x,2) + alpha * exact(x,1) + beta * exact(x,0);
}

} // namespace lawa
