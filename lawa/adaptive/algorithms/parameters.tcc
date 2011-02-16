namespace lawa {

template <typename T>
Parameters<T, Basis<T,Primal,R,CDF>, HelmholtzOperator1D<T,Basis<T,Primal,R,CDF> >,
           H1Preconditioner1D<T, Basis<T,Primal,R,CDF>, HelmholtzOperator1D<T,Basis<T,Primal,R,CDF> > > >::
           Parameters(const Basis<T,Primal,R,CDF> &_basis,
					  const HelmholtzOperator1D<T,Basis<T,Primal,R,CDF> > &_Bil,
					  bool _w_XBSpline, int _j0)
           : basis(_basis), Bil(_Bil), w_XBSpline(_w_XBSpline), j0(_j0),
             cA(0.), CA(0.), kappa(0.),
			 alpha(0.), omega(0.), gamma(0.), theta(0.)
{
    if (!_w_XBSpline) {	//only wavelet discretization
    	if (basis.d==2 && basis.d_==2 && Bil.c == 1.) {
    		cA = 0.19; CA = 2.;
    	}
    }
    else {				//wavelet with b-splines on coarse level discretization
    	if (j0==0) {
    		if (basis.d==2 && basis.d_==2 && Bil.c == 1.) {
    			cA = 0.375; CA = 2.1;

    		}
    	}
    }
    kappa = CA/cA;
	omega = 0.01;
	alpha = 1./std::sqrt(kappa)-(1.+1./std::sqrt(kappa))*omega-0.00001;
	gamma = 0.5 * (1./6.) * 1./sqrt(kappa) * (alpha-omega)/(1+omega);
	theta = 2./7.;

	assert(cA>0); assert(CA>0); assert(omega>0.); assert(alpha>0.); assert(gamma>0.); assert(theta>0.);
	assert((alpha+omega)/(1-omega)<1./std::sqrt(kappa));
}

template <typename T>
void
Parameters<T, Basis<T,Primal,R,CDF>, HelmholtzOperator1D<T,Basis<T,Primal,R,CDF> >,
           H1Preconditioner1D<T, Basis<T,Primal,R,CDF>, HelmholtzOperator1D<T,Basis<T,Primal,R,CDF> > > >::
           getGHSADWAVParameters(T &_alpha, T &_omega, T &_gamma, T &_theta) const
{
	_alpha = alpha; _omega = omega; _gamma = gamma; _theta = theta;
}

template <typename T>
void
Parameters<T, Basis<T,Primal,R,CDF>, HelmholtzOperator1D<T,Basis<T,Primal,R,CDF> >,
           H1Preconditioner1D<T, Basis<T,Primal,R,CDF>, HelmholtzOperator1D<T,Basis<T,Primal,R,CDF> > > >::
	getRHS_W_XBSplineParameters(T &_left_bound, T &_right_bound, int &_J_plus_smooth,
								 int &_J_plus_singular, bool &_singular_integral,
								 int example) const
{
	int d=basis.d, d_=basis.d_;
	if (example == 1) {
		_left_bound = -20.; _right_bound = 20.;
		_singular_integral=false;
		if (d==2 && d_==2) {
			_J_plus_smooth = 6;
			_J_plus_singular = 40;
		}
		else if (d==3) {
			_J_plus_smooth = 4;
			_J_plus_singular = 50;
		}
	}
	else if (example == 2) {
		_left_bound = -114.; _right_bound = 114.;
		_singular_integral=false;
		if (d==2 && d_==2) {
			_J_plus_smooth = 6;
			_J_plus_singular = 40;
		}
		else if (d==3) {
			_J_plus_smooth = 6;
			_J_plus_singular = 40;
		}
	}
	else if (example == 6) {
		_left_bound = -114.; _right_bound = 114.;
		_singular_integral=true;
		if (d==2 && d_==2) {
			_J_plus_smooth = 6;
			_J_plus_singular = 40;
		}
		else if (d==3) {
			_J_plus_smooth = 6;
			_J_plus_singular = 40;
		}
	}
	else {
		std::cerr << "Parameters not set for example " << example << " and d=" << d << ", d_=" << d_ << std::endl;
		assert(0);
	}

}


template <typename T>
void
Parameters<T, Basis<T,Primal,R,CDF>, HelmholtzOperator1D<T,Basis<T,Primal,R,CDF> >,
           H1Preconditioner1D<T, Basis<T,Primal,R,CDF>, HelmholtzOperator1D<T,Basis<T,Primal,R,CDF> > > >::
	getRHS_WO_XBSplineParameters(T &_left_bound, T &_right_bound, int &_J_plus_smooth, int &_J_minus_smooth,
								 int &_J_plus_singular, int &_J_minus_singular, bool &_singular_integral,
								 int example) const
{
	int d=basis.d, d_=basis.d_;
	if (example == 1) {
		_left_bound = -20.; _right_bound = 20.;
		_singular_integral=false;
		if (d==2 && d_==2) {
			_J_plus_smooth = 6; 	_J_minus_smooth = -40;
			_J_plus_singular = 40; 	_J_minus_singular = -40;
		}
		else if (d==3) {
			_J_plus_smooth = 4; 	_J_minus_smooth = -50;
			_J_plus_singular = 50;  _J_minus_singular = -50;
		}
	}
	else if (example == 2) {
		_left_bound = -114.; _right_bound = 114.;
		_singular_integral=false;
		if (d==2 && d_==2) {
			_J_plus_smooth = 6; 	_J_minus_smooth = -40;
			_J_plus_singular = 40; 	_J_minus_singular = -40;
		}
		else if (d==3) {
			_J_plus_smooth = 6; 	_J_minus_smooth = -40;
			_J_plus_singular = 40; 	_J_minus_singular = -40;
		}
	}
	else if (example == 6) {
		_left_bound = -114.; _right_bound = 114.;
		_singular_integral=true;
		if (d==2 && d_==2) {
			_J_plus_smooth = 6; 	_J_minus_smooth = -40;
			_J_plus_singular = 40; 	_J_minus_singular = -40;
		}
		else if (d==3) {
			_J_plus_smooth = 6; 	_J_minus_smooth = -40;
			_J_plus_singular = 40; 	_J_minus_singular = -40;
		}
	}
	else {
		std::cerr << "Parameters not set for example " << example << " and d=" << d << ", d_=" << d_ << std::endl;
		assert(0);
	}
}

}
