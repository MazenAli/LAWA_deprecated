namespace lawa {

template <typename T, typename Basis>
CGMYOperator1D<T, Basis>::CGMYOperator1D(const Basis& _basis, T _c, T _C, T _G, T _M, T _Y)
    : basis(_basis), c(_c), C(_C), G(_G), M(_M), Y(_Y),
      cgmy(C,G,M,Y),
      phi(basis.mra), d_phi(basis.mra, 1), psi(basis), d_psi(basis, 1),
      integral_sfsf(phi, phi), dd_integral_sfsf(d_phi, d_phi),
      integral_sfw(phi, psi), dd_integral_sfw(d_phi, d_psi),
      integral_wsf(psi, phi), dd_integral_wsf(d_psi, d_phi),
      integral_ww(psi,psi), dd_integral_ww(d_psi,d_psi)
{
	int jmin = basis.j0;
	int d    = basis.d;
	T step = 1.0/(1<<(jmin+7));
	Support<T> supp;
	phi_deltas.engine().resize(phi.singularSupport(jmin,0).length(),2);
	phi_deltas(_,1) = phi.singularSupport(jmin,0);
	supp = phi.support(jmin,0);

	PrimalSpline  dM1_th_deriv_phi(basis.mra,d-1);
	PrimalWavelet dM1_th_deriv_psi(basis,d-1);

	for (int i = 1; i<=phi_deltas.numRows(); ++i) {
		phi_deltas(i,2) = ((phi_deltas(i,1)==supp.l2) ? 0.0 : dM1_th_deriv_phi(std::min(phi_deltas(i,1)+step, supp.l2),jmin,0))
		        		- ((phi_deltas(i,1)==supp.l1) ? 0.0 : dM1_th_deriv_phi(std::max(phi_deltas(i,1)-step, supp.l1),jmin,0));

	}

	psi_deltas.engine().resize(psi.singularSupport(jmin,0).length(),2);
	psi_deltas(_,1) = psi.singularSupport(jmin,0);
	supp = psi.support(jmin,0);

    for (int i = 1; i<=psi_deltas.numRows(); ++i) {
    	psi_deltas(i,2) = ((psi_deltas(i,1)==supp.l2) ? 0.0 : dM1_th_deriv_psi(std::min(psi_deltas(i,1)+step, supp.l2),jmin,0))
    		    		- ((psi_deltas(i,1)==supp.l1) ? 0.0 : dM1_th_deriv_psi(std::max(psi_deltas(i,1)-step, supp.l1),jmin,0));
    }
}

template <typename T, typename Basis>
T
CGMYOperator1D<T,Basis>::getc() const
{
    return c;
}

template <typename T, typename Basis>
const Basis&
CGMYOperator1D<T,Basis>::getBasis() const
{
    return basis;
}

template <typename T, typename Basis>
GeMatrix<FullStorage<T,ColMajor> >
CGMYOperator1D<T, Basis>::computeDeltas(XType xtype, int j, int k) const
{
	int jmin=basis.j0;
	int d   =basis.d;
	GeMatrix<FullStorage<T,ColMajor> > ret;
	if (xtype == XBSpline) {
		ret.engine().resize(phi_deltas.numRows(),phi_deltas.numCols());
		for (int i = 1; i<=ret.numRows(); ++i) {
			ret(i,1) = phi_deltas(i,1)+pow2i<T>(-jmin)*k;
			ret(i,2) = phi_deltas(i,2);
		}
	}
	else {
		ret.engine().resize(psi_deltas.numRows(),psi_deltas.numCols());
		for (int i = 1; i<=ret.numRows(); ++i) {
			ret(i,1) = pow2i<T>(jmin-j)*psi_deltas(i,1)+pow2i<T>(-j)*k;
			ret(i,2) = pow2i<T>((d-1)*(j-jmin))*pow2ih<T>(j-jmin)*psi_deltas(i,2);
		}
	}
	return ret;
}

template <typename T, typename Basis>
T
CGMYOperator1D<T, Basis>::operator()(XType xtype1, int j1, int k1,
                                     XType xtype2, int j2, int k2) const
{
	//PDE part
    T val = 0;
    T dd_val = 0;

    if(xtype1 == XBSpline){
         if(xtype2 == XBSpline){
             val = integral_sfsf(j1, k1, j2, k2);
             dd_val = dd_integral_sfsf(j1, k1, j2, k2);
         }
         else{
             val = integral_sfw(j1, k1, j2, k2);
             dd_val = dd_integral_sfw(j1, k1, j2, k2);
         }
    }
    else{
         if(xtype2 == XBSpline){
             val = integral_wsf(j1, k1, j2, k2);
             dd_val = dd_integral_wsf(j1, k1, j2, k2);
         }
         else{
             val = integral_ww(j1, k1, j2, k2);
             dd_val = dd_integral_ww(j1, k1, j2, k2);
         }
    }

    GeMatrix<FullStorage<T,ColMajor> > psi_row_deltas = computeDeltas(xtype1,j1,k1);
    GeMatrix<FullStorage<T,ColMajor> > psi_col_deltas = computeDeltas(xtype2,j2,k2);
    T int_val = 0.;

    if (basis.d==2) {
    	for (int lambda=psi_row_deltas.rows().firstIndex(); lambda<=psi_row_deltas.rows().lastIndex(); ++lambda) {
    		T x = psi_row_deltas(lambda,1);
    		if (fabs(psi_row_deltas(lambda,2)) < 1e-14) continue;
    		if (xtype2 == XBSpline)	int_val += psi_row_deltas(lambda,2)*phi(x,j2,k2)*cgmy.c3;
    		else					int_val += psi_row_deltas(lambda,2)*psi(x,j2,k2)*cgmy.c3;

    		for (int mu=psi_col_deltas.rows().firstIndex(); mu<=psi_col_deltas.rows().lastIndex(); ++mu) {
    			if (fabs(psi_col_deltas(mu,2)) < 1e-14) continue;
    			T y = psi_col_deltas(mu,1);
    			T c = psi_col_deltas(mu,2)*psi_row_deltas(lambda,2);
    			if (x!=y)  {
    				//int_val += c * cgmy.nthTailIntegral(y-x,2*d,ZeroAtZero);
    				if (y-x>0)	int_val += c * (cgmy.ForthTailIntegral(y-x) - cgmy.constants[2]);
    				else 		int_val += c * (cgmy.ForthTailIntegral(y-x) - cgmy.constants[3]);
    			}
    		}
    	}
    }

	else if (basis.d==3) {

		if 		((xtype1 == XBSpline) && (xtype2 == XBSpline)) 	int_val -= cgmy.c3*dd_integral_sfsf(j1,k1,j2,k2);
		else if ((xtype1 == XBSpline) && (xtype2 == XWavelet))	int_val -= cgmy.c3*dd_integral_sfw(j1,k1,j2,k2);
		else if ((xtype1 == XWavelet) && (xtype2 == XBSpline))	int_val -= cgmy.c3*dd_integral_wsf(j1,k1,j2,k2);
		else													int_val -= cgmy.c3*dd_integral_ww(j1,k1,j2,k2);

		for (int mu=psi_col_deltas.rows().firstIndex(); mu<=psi_col_deltas.rows().lastIndex(); ++mu) {
			T y = psi_col_deltas(mu,1);
			if (fabs(psi_col_deltas(mu,2)) < 1e-14) continue;

			if (xtype1 == XBSpline)	{
				int_val += psi_col_deltas(mu,2)*phi(y,j1,k1)*cgmy.c4;
				int_val -= psi_col_deltas(mu,2)*d_phi(y,j1,k1)* cgmy.c5;
			}
			else {
				int_val += psi_col_deltas(mu,2)*psi(y,j1,k1)*cgmy.c4;
				int_val -= psi_col_deltas(mu,2)*d_psi(y,j1,k1)* cgmy.c5;
			}

			for (int lambda=psi_row_deltas.rows().firstIndex(); lambda<=psi_row_deltas.rows().lastIndex(); ++lambda) {
				if (fabs(psi_row_deltas(lambda,2)) < 1e-14) continue;
				T x = psi_row_deltas(lambda,1);
				T c = psi_col_deltas(mu,2)*psi_row_deltas(lambda,2);
				if (x!=y)  {
					if (y-x>0)		int_val -= c * (cgmy.SixthTailIntegral(y-x) - cgmy.constants[6]);
					else 			int_val -= c * (cgmy.SixthTailIntegral(y-x) - cgmy.constants[7]);
				}
			}
		}
	}
    else {
    	assert(0);
    }

    return dd_val +  c * val + int_val;
}

template <typename T, typename Basis>
T
CGMYOperator1D<T, Basis>::operator()(const Index1D &row_index, const Index1D &col_index) const
{
	return CGMYOperator1D<T, Basis>::operator()(row_index.xtype, row_index.j, row_index.k,
											    col_index.xtype, col_index.j, col_index.k);
}

}	//namespace lawa
