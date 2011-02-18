namespace lawa {

template <typename T, typename Index, typename Basis1D, typename APPLY1D, typename RHS>
GHS_ADWAV1D<T,Index,Basis1D,APPLY1D,RHS>::GHS_ADWAV1D(const Basis1D &_basis, APPLY1D &_Apply,
													  RHS &_F)
	: basis(_basis), Apply(_Apply), F(_F),
	  cA(0.), CA(0.), kappa(0.),
	  alpha(0.), omega(0.), gamma(0.), theta(0.), eps(0.)
{
	cA    = Apply.parameters.cA;
	CA    = Apply.parameters.CA;
	kappa = Apply.parameters.kappa;
	Apply.parameters.getGHSADWAVParameters(alpha, omega, gamma, theta);


}

template <typename T, typename Index, typename Basis1D, typename APPLY1D, typename RHS>
Coefficients<Lexicographical,T,Index>
GHS_ADWAV1D<T,Index,Basis1D,APPLY1D,RHS>::SOLVE(T nuM1, T _eps, int NumOfIterations)
{
	T eps = _eps;
	int d=basis.d, d_=basis.d_;
	Coefficients<Lexicographical,T,Index> w_k(d,d_), w_kP1(d,d_), g_kP1(d,d_);
	IndexSet<Index> Lambda_kP1(d,d_);
	T nu_kM1 = nuM1;
	T nu_k   = 0.;
	T total_time=0.;

	std::cout << "SOLVE has started with the following parameters: " << std::endl;
	std::cout << "  alpha=" << alpha << ", omega=" << omega << ", gamma=" << gamma << ", theta="
			  << theta << std::endl;
	std::cout << "  cA=" << cA << ", CA=" << CA << ", kappa=" << kappa << std::endl;

	for (int i=1; i<=NumOfIterations; ++i) {
		Timer time;
		std::cout << "*** " << i << ".iteration ***" << std::endl;
		time.start();
		std::cout << "  GROW started." << std::endl;
		Lambda_kP1 = GROW(w_k, theta*nu_kM1, nu_k);
		std::cout << "  GROW finished." << std::endl;
		solutions.push_back(w_kP1);
		residuals.push_back(nu_k);
		times.push_back(total_time);
		if (nu_k <=eps) break;

		std::cout << "   GALSOLVE started with #Lambda = " << Lambda_kP1.size()  << std::endl;
		//g_kP1 = P(F(gamma*nu_k),Lambda_kP1);
		g_kP1 = F(Lambda_kP1);
		w_kP1 = GALSOLVE(Lambda_kP1, g_kP1, w_k, (1+gamma)*nu_k, gamma*nu_k);
		//T r_norm_LambdaActive;
		//int iterations = CG_Solve(Lambda_kP1, Apply.A, w_kP1, g_kP1, r_norm_LambdaActive, 1e-16);
		//std::cout << "   iterations = " << iterations << ", residual = " << r_norm_LambdaActive << ", w_k = " << w_k << std::endl;
		std::cout << "  GALSOLVE finished." << std::endl;
		nu_kM1 = nu_k;
		w_k = w_kP1;
		time.stop();
		total_time += time.elapsed();
		std::cout << std::endl;
	}
	return w_k;
}

template <typename T, typename Index, typename Basis1D, typename APPLY1D, typename RHS>
IndexSet<Index>
GHS_ADWAV1D<T,Index,Basis1D,APPLY1D,RHS>::GROW(const Coefficients<Lexicographical,T,Index> &w, T nu_bar, T &nu)
{
	int d=basis.d, d_=basis.d_;
	T zeta = 2.*(omega*nu_bar)/(1-omega);
	T r_norm = 0.;
	Coefficients<Lexicographical,T,Index> r(d,d_), Aw(d,d_), rhs(d,d_);
	while (1) {
		zeta *= 0.5;
		rhs = F(0.5*zeta);
		Aw  = Apply(w, 0.5*zeta);
		r = rhs - Aw;
		r_norm = r.norm(2.);
		nu = r_norm + zeta;
		//std::cout << "    zeta = " << zeta << ", r_norm = " << r_norm << ", omega*r_norm = " << omega*r_norm << ", nu = " << nu << std::endl;
		if (nu <= eps) break;
		if (zeta<=omega*r_norm) break;
	}


	Coefficients<AbsoluteValue,T,Index> r_abs(d,d_);
	r_abs = r;
	IndexSet<Index1D> Lambda(d,d_);
	T P_Lambda_r_norm_square = 0.;
	if (w.size() > 0) {
		Lambda = supp(w);
		for (const_coeff_it it=w.begin(); it!=w.end(); ++it) {
			P_Lambda_r_norm_square += std::pow(r[(*it).first],2.);
		}
	}

	//std::cout << "   Before extension: ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_r_norm_square) << ", alpha*r_norm = " << alpha*r_norm << std::endl;
	if (nu > eps) {
		for (const_coeff_abs_it it=r_abs.begin(); it!=r_abs.end(); ++it) {
			if (Lambda.count((*it).second) == 0) {
				Lambda.insert((*it).second);
				P_Lambda_r_norm_square += std::pow((*it).first,2);
				//std::cout << "    Added " << (*it).second << ", now: ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_r_norm_square) << ", alpha*r_norm = " << alpha*r_norm << std::endl;
				if (P_Lambda_r_norm_square >= alpha*r_norm*alpha*r_norm) break;
			}
		}
	}
	return Lambda;
}

template <typename T, typename Index, typename Basis1D, typename APPLY1D, typename RHS>
Coefficients<Lexicographical,T,Index>
GHS_ADWAV1D<T,Index,Basis1D,APPLY1D,RHS>::GALSOLVE(const IndexSet<Index> &Lambda,
											       const Coefficients<Lexicographical,T,Index> &g,
											       const Coefficients<Lexicographical,T,Index> &w,
											       T delta, T tol)
{
	int d=basis.d, d_=basis.d_;
	Coefficients<Lexicographical,T,Index> ret(d,d_);

	if (Lambda.size()==0) return ret;

	//Determine compression level
	int J=0;		//compression
	if 		(d==2) {   J = -std::ceil(2*std::log(tol/((3*tol+3*delta)*kappa))); }
	else if (d==3) {   J = -std::ceil((1./1.5)*std::log(tol/((3*tol+3*delta)*kappa))); }
	else 			 { assert(0); }
	//std::cout << "    Estimated compression level for delta=" << delta << " and target tol=" << tol
	//		  << " : " << J << std::endl;

	//Assemble sparse matrix B
	int N = Lambda.size();
	//std::cout << "    Assembling of B started with N=" << N << std::endl;
	flens::SparseGeMatrix<CRS<T,CRS_General> > B(N,N);
	std::map<Index,int,lt<Lexicographical,Index> > row_indices;
    int row_count = 1, col_count = 1;
    for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
    	row_indices[(*row)] = row_count;
    }
    Apply.A.c.setParameters(Lambda);
    for (const_set_it col=Lambda.begin(); col!=Lambda.end(); ++col, ++col_count) {
    	IndexSet<Index> LambdaRowSparse = Apply.A.c.SparsityPattern(*col, Lambda, J);
    	for (const_set_it row=LambdaRowSparse.begin(); row!=LambdaRowSparse.end(); ++row) {
    		T tmp = Apply.A(*row,*col);
    		if (fabs(tmp)>0)	B(row_indices[*row],col_count) = 0.5*tmp;
    		if (fabs(tmp)>0)	B(col_count,row_indices[*row]) = 0.5*tmp;
    	}
    }
    B.finalize();
    //std::cout << "    Assembling of B finished." << std::endl;

    //Compute approximate initial residual
    //std::cout << "    Solving Bx=r0." << std::endl;
	Coefficients<Lexicographical,T,Index> r0(d,d_), APPLY_Aw(d,d_);
	APPLY_Aw = Apply(w, tol/3.);
	r0 = g - P(APPLY_Aw, Lambda);

	DenseVector<Array<T> > rhs(N), x(N), res(N), Bx(N);
	row_count=1;
	const_coeff_it r0_end = r0.end();
	for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
		const_coeff_it it = r0.find(*row);
		if (it != r0_end) rhs(row_count) = (*it).second;
	    else                      rhs(row_count) = 0.;
	}
	//std::cout << "    cg-method started with rhs " << rhs << std::endl;
	int iters = lawa::cg(B,x,rhs,tol);
	Bx = B*x;
	res= Bx-rhs;
	T lin_res = std::sqrt(res*res);
	std::cout << "    cg-method needed " << iters << " iterations, residual=" << lin_res << std::endl;
	assert(lin_res<tol);
	row_count = 1;
	const_coeff_it w_end = w.end();
	for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
		T tmp = 0.;
		const_coeff_it it = w.find(*row);
		if (it != w_end) tmp = (*it).second;
		ret[*row] = tmp+x(row_count);
	}
	//std::cout << "    Calculated Bx=r0, w+x=" << ret << std::endl;
	return ret;
}


}	//namespace lawa

