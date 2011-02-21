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
    	 std::cout << "Parameters: j0 = - Infinity" << std::endl;
    	if (basis.d==2 && basis.d_==2 && Bil.c == 1.) {
    		cA = 0.19; CA = 2.;
    	}
    	else if (basis.d==3 && basis.d_==3 && Bil.c == 1.) {
    		cA = 0.038; CA = 2.7;
    	}
    	else if (basis.d==3 && basis.d_==5 && Bil.c == 1.) {
    		cA = 0.16; CA = 2.8;
    	}
    	else assert(0);
    }
    else {				//wavelet with b-splines on coarse level discretization
    	 std::cout << "Parameters: j0 = " << j0 << std::endl;
    	if (basis.d==2 && basis.d_==2 && Bil.c == 1.) {
    		if (j0==0)  {	cA = 0.375; CA = 2.1;	}
    		else if (j0==-1) {	cA = 0.58;  CA = 1.86;	}
    		else if (j0==-2) {	cA = 0.58;  CA = 1.86;	}
    		else if (j0==-3) {	cA = 0.55;  CA = 1.86;	}
    		else if (j0==-4) {	cA = 0.46;  CA = 1.86;	}
    		else if (j0==-5) {	cA = 0.4;   CA = 1.86;	}
    		else if (j0==-6) {	cA = 0.36;  CA = 1.89;	}
    		else if (j0==-7) {	cA = 0.33;  CA = 1.92;	}
    		else if (j0==-8) {	cA = 0.3;   CA = 1.95;	}
    		else if (j0==-8) {	cA = 0.29;  CA = 1.95;	}
    		else if (j0==-10) {	cA = 0.27;  CA = 1.96;	}
    		else assert(0);
    	}
    	else if (basis.d==3 && basis.d_==3 && Bil.c == 1.) {
    	    if (j0==0)  {	cA = 0.43;  CA = 1.94;	}
    	    else if (j0==-1) {	cA = 0.39;  CA = 2.03;	}
    	    else if (j0==-2) {	cA = 0.29;  CA = 2.4;	}
    	    else if (j0==-3) {	cA = 0.21;  CA = 2.44;	}
    	    else if (j0==-4) {	cA = 0.16;  CA = 2.55;	}
    	    else if (j0==-5) {	cA = 0.12;  CA = 2.59;	}
    	    else if (j0==-6) {	cA = 0.1;   CA = 2.61;	}
    	    else if (j0==-7) {	cA = 0.09;  CA = 2.62;	}
    	    else if (j0==-8) {	cA = 0.078; CA = 2.63;	}
    	    else if (j0==-9) {	cA = 0.07;  CA = 2.64;	}
    	    else if (j0==-10) {	cA = 0.063; CA = 2.64;	}
    	    else assert(0);
    	}
    	else if (basis.d==3 && basis.d_==5 && Bil.c == 1.) {
    		if (j0==0)  {	cA = 0.457;  CA = 1.95;	}
    		else if (j0==-1) {	cA = 0.416;  CA = 2.07;	}
    		else if (j0==-2) {	cA = 0.326;  CA = 2.32;	}
    	    else if (j0==-3) {	cA = 0.24;   CA = 2.53;	}
    	    else if (j0==-4) {	cA = 0.197;  CA = 2.66;	}
    	    else if (j0==-5) {	cA = 0.18;   CA = 2.7;	}
    	    else assert(0);
    	}
    	else assert(0);
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

		_singular_integral=false;
		if (d==2 && d_==2) {
			_left_bound = -20.; _right_bound = 20.;
			_J_plus_smooth = 6;
			_J_plus_singular = 40;
		}
		else if (d==3) {
			_left_bound = -80.; _right_bound = 80.;
			_J_plus_smooth = 4;
			_J_plus_singular = 50;
		}
	}
	else if (example == 2) {
		_singular_integral=true;
		if (d==2 && d_==2) {
			_left_bound = -114.; _right_bound = 114.;
			_J_plus_smooth = 6;
			_J_plus_singular = 40;
		}
		else if (d==3) {
			_left_bound = -214.; _right_bound = 214.;
			_J_plus_smooth = 5;
			_J_plus_singular = 40;
		}
	}
	else if (example == 3) {
		_singular_integral=true;
		if (d==2 && d_==2) {
			_left_bound = -100.; _right_bound = 100.;
			_J_plus_smooth = 6;
			_J_plus_singular = 40;
		}
		else if (d==3) {
			_left_bound = -100.; _right_bound = 100.;
			_J_plus_smooth = 4;
			_J_plus_singular = 60;
		}
	}
	else if (example == 4) {
		_singular_integral=true;
		_left_bound = -0.4; _right_bound = 0.9;
		if (d==2 && d_==2) {
			_J_plus_smooth = 6;
			_J_plus_singular = 40;
		}
		else if (d==3) {
			_J_plus_smooth = 5;
			_J_plus_singular = 50;
		}
	}
	else if (example == 5) {
		_singular_integral=true;
		_left_bound = -0.4; _right_bound = 0.9;
		if (d==2 && d_==2) {
			_J_plus_smooth = 8;
			_J_plus_singular = 40;
		}
		else if (d==3) {
			_J_plus_smooth = 6;
			_J_plus_singular = 40;
		}
	}
	else if (example == 6) {
		_singular_integral=true;
		if (d==2 && d_==2) {
			_left_bound = -114.; _right_bound = 114.;
			_J_plus_smooth = 6;
			_J_plus_singular = 40;
		}
		else if (d==3) {
			_left_bound = -200.; _right_bound = 200.;
			_J_plus_smooth = 3;
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
		_singular_integral=false;
		if (d==2 && d_==2) {
			_left_bound = -20.; _right_bound = 20.;
			_J_plus_smooth = 6; 	_J_minus_smooth = -40;
			_J_plus_singular = 40; 	_J_minus_singular = -40;
		}
		else if (d==3) {
			_left_bound = -40.; _right_bound = 40.;
			_J_plus_smooth = 4; 	_J_minus_smooth = -50;
			_J_plus_singular = 50;  _J_minus_singular = -50;
		}
	}
	else if (example == 2) {
		_singular_integral=true;
		if (d==2 && d_==2) {
			_left_bound = -114.; _right_bound = 114.;
			_J_plus_smooth = 6; 	_J_minus_smooth = -40;
			_J_plus_singular = 40; 	_J_minus_singular = -40;
		}
		else if (d==3) {
			_left_bound = -200.; _right_bound = 200.;
			_J_plus_smooth = 4; 	_J_minus_smooth = -80;
			_J_plus_singular = 80; 	_J_minus_singular = -80;
		}
	}
	else if (example == 3) {
		_singular_integral=true;
		if (d==2 && d_==2) {
			_left_bound = -50.; _right_bound = 50.;
			_J_plus_smooth = 6; 	_J_minus_smooth = -40;
			_J_plus_singular = 40; 	_J_minus_singular = -40;
		}
		else if (d==3) {
			_left_bound = -100.; _right_bound = 100.;
			_J_plus_smooth = 3; 	_J_minus_smooth = -80;
			_J_plus_singular = 80; 	_J_minus_singular = -80;
		}
	}
	else if (example == 4) {
		_singular_integral=true;
		_left_bound = -0.4; _right_bound = 0.9;
		if (d==2 && d_==2) {
			_J_plus_smooth = 6; 	_J_minus_smooth = -40;
			_J_plus_singular = 40; 	_J_minus_singular = -40;
		}
		else if (d==3) {
			_J_plus_smooth = 5; 	_J_minus_smooth = -50;
			_J_plus_singular = 50; 	_J_minus_singular = -50;
		}
	}
	else if (example == 5) {
		_singular_integral=true;
		_left_bound = -0.4; _right_bound = 0.9;
		if (d==2 && d_==2) {
			_J_plus_smooth = 10; 	_J_minus_smooth = -40;
			_J_plus_singular = 40; 	_J_minus_singular = -40;
		}
		else if (d==3) {
			_J_plus_smooth = 9; 	_J_minus_smooth = -50;
			_J_plus_singular = 45; 	_J_minus_singular = -50;
		}
	}
	else if (example == 6) {
		_singular_integral=true;
		if (d==2 && d_==2) {
			_left_bound = -114.; _right_bound = 114.;
			_J_plus_smooth = 6; 	_J_minus_smooth = -40;
			_J_plus_singular = 40; 	_J_minus_singular = -40;
		}
		else if (d==3) {
			_left_bound = -200.; _right_bound = 200.;
			_J_plus_smooth = 3; 	_J_minus_smooth = -40;
			_J_plus_singular = 40; 	_J_minus_singular = -40;
		}
	}
	else {
		std::cerr << "Parameters not set for example " << example << " and d=" << d << ", d_=" << d_ << std::endl;
		assert(0);
	}
}

}
