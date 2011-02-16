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
RHSWithPeaks1D<T,Basis>::RHSWithPeaks1D(bool _with_singular_part, bool _with_smooth_part,
										const Basis &_basis, T (*_f)(T),
										const DenseVector<Array<T> > &_singularPoints,
										const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &_deltas,
										int order)
    : with_singular_part(_with_singular_part), with_smooth_part(_with_smooth_part),
	  basis(_basis),  phi(basis.mra), psi(basis),
      f(_f,_singularPoints), deltas(_deltas),
      integral_sff(phi, f), integral_wf(psi, f)
{
    integral_sff.quadrature.setOrder(order);
    integral_wf.quadrature.setOrder(order);
}

template <typename T, typename Basis>
T
RHSWithPeaks1D<T,Basis>::operator()(XType xtype, int j, int k) const
{
    T ret = 0.;
    if (xtype == XBSpline)  {
        if (with_smooth_part) ret = integral_sff(j,k);
        if (with_singular_part) {
        	for (int i=1; i<=deltas.numRows(); ++i) {
        		ret += deltas(i,2) * phi(deltas(i,1),j,k);
        	}
        }
    }
    else {
    	if (with_smooth_part) ret = integral_wf(j,k);
    	if (with_singular_part) {
    		for (int i=1; i<=deltas.numRows(); ++i) {
    			ret += deltas(i,2) * psi(deltas(i,1),j,k);
    		}
    	}
    }

    return ret;
}

template <typename T, typename Basis>
T
RHSWithPeaks1D<T,Basis>::operator()(const Index1D &lambda) const
{
    return RHSWithPeaks1D<T,Basis>::operator()(lambda.xtype, lambda.j, lambda.k);
}



template <typename T>
RHSWithPeaks1D_WO_XBSpline<T>::RHSWithPeaks1D_WO_XBSpline(bool _with_singular_part, bool _with_smooth_part,
														  const Wavelet<T,Primal,R,CDF> &_psi,
														  T _left_bound, T _right_bound,  T (*_f)(T),
														  const DenseVector<Array<T> > &_f_singularPoints,
														  flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > &_deltas,
														  T _h, int _order)
: with_singular_part(_with_singular_part), with_smooth_part(_with_smooth_part), psi(_psi),
  f(_f), left_bound(_left_bound), right_bound(_right_bound), f_singularPoints(_f_singularPoints),
  deltas(_deltas), h(_h), order(_order)
{
	DenseVector<Array<T> > bounds(2);
	bounds(1) = left_bound;
	bounds(2) = right_bound;

	int nFirst  = bounds.length(),
		nSecond = f_singularPoints.length();

	f_singularPoints_interval.engine().resize(nFirst + nSecond);

	std::merge(bounds.engine().data(),
		       bounds.engine().data() + nFirst,
		       f_singularPoints.engine().data(),
		       f_singularPoints.engine().data() + nSecond,
		       f_singularPoints_interval.engine().data());


	T eps = Const<T>::EQUALITY_EPS;
    _knots.engine().resize(_order, _order);
    _weights.engine().resize(_order, _order);

    T x1 = -1,
      x2 =  1;
    for (int k=1; k<=_order; ++k) {
        int     m = (k+1)/2;
        T xm = 0.5 * (x2+x1),
          xl = 0.5 * (x2-x1);
        for (int i=1; i<=m; ++i) {
            T z = cos(M_PI*(i-0.25)/(k+0.5)),
              z1, pp;
            do {
                T p1 = 1.0,
                  p2 = 2.0;
                for (int j=1; j<=k; ++j) {
                    T p3 = p2;
                    p2 = p1;
                    p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
                }
                pp = k * (z*p1-p2)/(z*z-1.0);
                z1 = z;
                z = z1-p1/pp;
            } while (fabs(z-z1) > eps);
            _knots(k,i)     = xm - xl*z;
            _knots(k,k+1-i) = xm + xl*z;
            _weights(k,i)     = 2.0*xl/((1.0-z*z)*pp*pp);
            _weights(k,k+1-i) = _weights(k,i);
        }
    }
}

template <typename T>
T
RHSWithPeaks1D_WO_XBSpline<T>::truncated_f(T x) const
{
	if 		(x <= left_bound-1) return 0.;
	else if ((left_bound-1 < x) && (x<left_bound)) {
		T help1 = (-left_bound+1+x)*(-left_bound+1+x);
		T help2 = (-left_bound+x)*(-left_bound+x);
		return f(x)*std::exp(-1./help1)/(std::exp(-1./help1)+std::exp(-1./help2));
	}
	else if ((left_bound <= x) && (x<=right_bound)) return f(x);
	else if ((right_bound < x) && (x<right_bound+1)) {
		T help1 = (right_bound+1-x)*(right_bound+1-x);
		T help2 = (right_bound-x)*(right_bound-x);
		return f(x)*std::exp(-1./help1)/(std::exp(-1./help1)+std::exp(-1./help2));
	}
	else return 0.;
}

template <typename T>
void
RHSWithPeaks1D_WO_XBSpline<T>::computeWaveletNorms() const
{
	T psi_L1norm = 0., psi_L2norm = 0., d_psi_L2norm = 0.;
	int d = psi.d;
	int d_= psi.d_;
	int int_order = 2*d;
	DenseVector<Array<T> > psi_singularPoints = psi.singularSupport(0,0);
	for (int j=psi_singularPoints.firstIndex(); j<psi_singularPoints.lastIndex(); ++j) {
		T a = psi_singularPoints(j), b = psi_singularPoints(j+1);
		T tmp1 = 0., tmp2 = 0.;
		for (int i=1; i<=order; ++i) {
		    T x = 0.5*(b-a)*_knots(int_order,i)+0.5*(b+a);
		    tmp1 += _weights(int_order,i) * fabs(psi(x,0,0));
		    tmp2 += _weights(int_order,i) * psi(x,0,0) * psi(x,0,0);
		}
		psi_L1norm += 0.5*(b-a)*tmp1;
		psi_L2norm += 0.5*(b-a)*tmp2;
	}
	if (d>1) {
		Wavelet<T,Primal,R,CDF> d_psi(d-1,d_+1);
		DenseVector<Array<T> > d_psi_singularPoints = d_psi.singularSupport(0,0);
		for (int j=d_psi_singularPoints.firstIndex(); j<d_psi_singularPoints.lastIndex(); ++j) {
			T a = d_psi_singularPoints(j), b = d_psi_singularPoints(j+1);
			T tmp = 0.;
			for (int i=1; i<=order; ++i) {
				T x = 0.5*(b-a)*_knots(int_order,i)+0.5*(b+a);
				tmp += _weights(int_order,i) * d_psi(x,0,0) * d_psi(x,0,0);
			}
			d_psi_L2norm += 0.5*(b-a)*tmp;
		}
		std::cout << "||psi||_L1   = " << psi_L1norm << std::endl;
		std::cout << "||psi||_L2   = " << psi_L2norm << std::endl;
		std::cout << "||d_psi||_L2 = " << d_psi_L2norm << std::endl;
	}
	else {
		std::cout << "||psi||_L1   = " << psi_L1norm << std::endl;
		std::cout << "||psi||_L2   = " << psi_L2norm << std::endl;
	}
}

template <typename T>
T
RHSWithPeaks1D_WO_XBSpline<T>::operator()(const Index1D &lambda) const
{
	//std::cout << "Integral_WO_XBSpline<T>::operator(" << lambda << ") was called." << std::endl;

	T ret = 0.0;

	if (with_singular_part) {
		for (int i=1; i<=deltas.numRows(); ++i) {
			ret += deltas(i,2) * psi(deltas(i,1),lambda.j,lambda.k);
		}
	}

	if (with_smooth_part) {
		DenseVector<Array<T> > psi_singularPoints = psi.optim_singularSupport(lambda.j,lambda.k);

		if (lambda.j >= 0) {
			int nFirst  = psi_singularPoints.length(),
				nSecond = f_singularPoints.length();

			DenseVector<Array<T> > singularPoints(nFirst + nSecond);

			std::merge(psi_singularPoints.engine().data(),
					   psi_singularPoints.engine().data() + nFirst,
					   f_singularPoints.engine().data(),
					   f_singularPoints.engine().data() + nSecond,
					   singularPoints.engine().data());


			for (int i=singularPoints.firstIndex(); i<singularPoints.lastIndex(); ++i) {
				ret += RHSWithPeaks1D_WO_XBSpline<T>::operator()(lambda.j,lambda.k,singularPoints(i),singularPoints(i+1));
			}
		}
		else {
			int nFirst  = psi_singularPoints.length(),
				nSecond = f_singularPoints_interval.length();

			DenseVector<Array<T> > singularPoints(nFirst + nSecond);

			std::merge(psi_singularPoints.engine().data(),
					   psi_singularPoints.engine().data() + nFirst,
					   f_singularPoints_interval.engine().data(),
					   f_singularPoints_interval.engine().data() + nSecond,
					   singularPoints.engine().data());

			T left_integral_bound  = std::max(left_bound-1, psi.support(lambda.j,lambda.k).l1);
			T right_integral_bound = std::min(right_bound+1,psi.support(lambda.j,lambda.k).l2);

			for (int i=singularPoints.firstIndex(); i<singularPoints.lastIndex(); ++i) {

				if      (singularPoints(i+1)<=left_integral_bound)  continue;
				else if (singularPoints(i)  >=right_integral_bound) break;
				else {
					T left, right;
					if (singularPoints(i+1) >left_integral_bound && singularPoints(i)<=left_integral_bound) {
						left=left_integral_bound;
						right=singularPoints(i+1);
					}
					else if (singularPoints(i+1)<=right_integral_bound && singularPoints(i)>=left_integral_bound) {
						left=singularPoints(i);
						right=singularPoints(i+1);
					}
					else if (singularPoints(i+1)> right_integral_bound && singularPoints(i)<=right_integral_bound) {
						left=singularPoints(i);
						right=right_integral_bound;
					}
					else {
						std::cerr << "Integral_WO_XBSpline<T>::operator()" << std::endl;
						exit(1);
					}
					T stepsize = min(h,right-left);
					T a = left;
					while (a<right) {
						T b = a+stepsize;
						//std::cout << "[" << a << ", " << b << "]" << std::endl;
						ret += RHSWithPeaks1D_WO_XBSpline<T>::operator()(lambda.j,lambda.k,a,b);
						a = b;
						stepsize = min(h,right-b);
					}
					if (fabs(a-right)>1e-14) {
						//std::cout << "rest: [" << a << ", " << right << "]" << std::endl;
						ret += RHSWithPeaks1D_WO_XBSpline<T>::operator()(lambda.j,lambda.k,a,right);
					}
				}
			}
		}
	}
	return ret;
}

template <typename T>
T
RHSWithPeaks1D_WO_XBSpline<T>::operator()(int j, int k, T a, T b) const
{
    T ret = 0.0;
    for (int i=1; i<=order; ++i) {
    	T x = 0.5*(b-a)*_knots(order,i)+0.5*(b+a);
        ret += _weights(order,i) * psi(x,j,k) * truncated_f(x);
    }
    ret *= 0.5*(b-a);

    return ret;
}



}  //namespace lawa
