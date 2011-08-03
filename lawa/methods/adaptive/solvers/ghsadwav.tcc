namespace lawa {

template <typename T, typename Index, typename AdaptiveOperator, typename RHS>
GHS_ADWAV<T,Index,AdaptiveOperator,RHS>::GHS_ADWAV(AdaptiveOperator &_A, RHS &_F)
    : A(_A), F(_F),
      cA(A.cA), CA(A.CA), kappa(A.kappa),
      alpha(0.), omega(0.), gamma(0.), theta(0.), eps(0.),
      sparseMat_A(200000,200000)
{
    omega = 0.01;
    alpha = 1./std::sqrt(kappa)-(1.+1./std::sqrt(kappa))*omega-0.00001;
    gamma = 0.5 * (1./6.) * 1./sqrt(kappa) * (alpha-omega)/(1+omega);
    theta = 2./7.;

}

template <typename T, typename Index, typename AdaptiveOperator, typename RHS>
Coefficients<Lexicographical,T,Index>
GHS_ADWAV<T,Index,AdaptiveOperator,RHS>::SOLVE(T nuM1, T _eps, int NumOfIterations, T H1norm)
{
    T eps = _eps;
    Coefficients<Lexicographical,T,Index> w_k, w_kP1, g_kP1;
    IndexSet<Index> Lambda_kP1;
    T nu_kM1 = nuM1;
    T nu_k   = 0.;
    T total_time=0.;

    std::cerr << "GHS-ADWAV-SOLVE has started with the following parameters: " << std::endl;
    std::cerr << "  alpha=" << alpha << ", omega=" << omega << ", gamma=" << gamma << ", theta="
              << theta << std::endl;
    std::cerr << "  cA=" << cA << ", CA=" << CA << ", kappa=" << kappa << std::endl;

    std::stringstream filename;
    filename << "ghs-adwav-nd-otf.dat";
    std::ofstream file(filename.str().c_str());

    for (int i=1; i<=NumOfIterations; ++i) {
        Timer time;
        std::cerr << "*** " << i << ".iteration ***" << std::endl;
        time.start();
        std::cerr << "  GROW started." << std::endl;

        IndexSet<Index> Extension;
        Extension = this->GROW(w_k, theta*nu_kM1, nu_k);
        Lambda_kP1 = Lambda_kP1 + Extension;
        std::cerr << "  GROW finished." << std::endl;
        //solutions.push_back(w_kP1);
        residuals.push_back(nu_k);
        times.push_back(total_time);
        if (nu_k <=eps) break;

        time.stop();
        total_time += time.elapsed();

        Coefficients<Lexicographical,T,Index> Au, f;
        f = F(supp(w_k));
        T fu = w_k*f;
        Au = A.mv(supp(w_k), w_k);
        T uAu = w_k*Au;
        T Error_H_energy = sqrt(fabs(std::pow(H1norm,2.)- 2*fu + uAu));
        //T Error_H_energy = computeErrorInH1Norm(Apply.A, F, w_k, H1norm);
        file << w_k.size() << " " << total_time << " " <<  nu_k << " "
                         << Error_H_energy << std::endl;


        time.start();
        std::cerr << "   GALSOLVE started with #Lambda = " << Lambda_kP1.size()  << std::endl;
        //g_kP1 = F(Lambda_kP1);                // update rhs vector
        g_kP1 = P(F(gamma*nu_k),Lambda_kP1);  // compute with restriction, otherwise GROW does not work
        w_kP1 = this->GALSOLVE(Lambda_kP1, Extension, g_kP1, w_k, (1+gamma)*nu_k, gamma*nu_k);
        //T r_norm_LambdaActive;
        //int iterations = CG_Solve(Lambda_kP1, Apply.A, w_kP1, g_kP1, r_norm_LambdaActive, 1e-16);
        //std::cerr << "   iterations = " << iterations << ", residual = "
        //          << r_norm_LambdaActive << ", w_k = " << w_k << std::endl;
        std::cerr << "  GALSOLVE finished." << std::endl;

        nu_kM1 = nu_k;
        w_k = w_kP1;
        time.stop();
        total_time += time.elapsed();
        std::cerr << std::endl;
    }
    return w_k;
}

template <typename T, typename Index, typename AdaptiveOperator, typename RHS>
IndexSet<Index>
GHS_ADWAV<T,Index,AdaptiveOperator,RHS>::GROW(const Coefficients<Lexicographical,T,Index> &w,
                                              T nu_bar, T &nu)
{
    T zeta = 2.*(omega*nu_bar)/(1-omega);
    T r_norm = 0.;
    Coefficients<Lexicographical,T,Index> r, Aw, rhs;
    while (1) {
        zeta *= 0.5;
        rhs = F(0.5*zeta);
        Aw  = A.apply(w, 0.5*zeta);
        r = rhs - Aw;
        r_norm = r.norm(2.);
        nu = r_norm + zeta;
        std::cerr << "    zeta = " << zeta << ", r_norm = " << r_norm
                  << ", omega*r_norm = " << omega*r_norm << ", nu = " << nu << std::endl;
        if (nu <= eps) break;
        if (zeta<=omega*r_norm) break;
    }

    IndexSet<Index> Lambda;
    T P_Lambda_r_norm_square = 0.;
    if (w.size() > 0) {
        for (const_coeff_it it=w.begin(); it!=w.end(); ++it) {
            //Lambda.insert((*it).first);
            P_Lambda_r_norm_square += std::pow(r[(*it).first],2.);
            r.erase((*it).first);
        }
    }

    T threshbound = std::sqrt(1-alpha*alpha) * r.norm(2.)/std::sqrt(T(r.size()));
    Coefficients<Bucket,T,Index> r_bucket;
    r_bucket.bucketsort(r, threshbound);
    std::cerr << "   Threshbound = " << threshbound << std::endl;
    std::cerr << "   Before extension: ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_r_norm_square)
              << ", alpha*r_norm = " << alpha*r_norm << std::endl;
    if (nu > eps) {
        int currentDoF = w.size();
        for (int i=0; i<(int)r_bucket.bucket_ell2norms.size(); ++i) {
            P_Lambda_r_norm_square += std::pow(r_bucket.bucket_ell2norms[i],2.);
            std::cerr << "   ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_r_norm_square) << std::endl;
            int addDoF = r_bucket.addBucketToIndexSet(Lambda,i,currentDoF);
            currentDoF += addDoF;
            if (P_Lambda_r_norm_square >= alpha*r_norm*alpha*r_norm) break;
        }
        /*
        int count=0;
        int sizeExtensionOfLambda=1;
        for (const_coeff_abs_it it=r_abs.begin(); it!=r_abs.end(); ++it) {

            Lambda.insert((*it).second);
            ++count;

            P_Lambda_r_norm_square += std::pow((*it).first,2);
            std::cerr << "    Added " << (*it).second << ", " << (*it).first << ", now: ||P_{Lambda}r ||_2 = "
                      << std::sqrt(P_Lambda_r_norm_square) << ", alpha*r_norm = "
                      << alpha*r_norm << std::endl;

            if (count>=5*sizeExtensionOfLambda)                     break;
        }
        */
    }
    //getchar();
    return Lambda;
}

template <typename T, typename Index, typename AdaptiveOperator, typename RHS>
Coefficients<Lexicographical,T,Index>
GHS_ADWAV<T,Index,AdaptiveOperator,RHS>::GALSOLVE(const IndexSet<Index> &Lambda,
                                                  const IndexSet<Index> &Extension,
                                                  const Coefficients<Lexicographical,T,Index> &g,
                                                  const Coefficients<Lexicographical,T,Index> &w,
                                                  T delta, T tol)
{
    int d=A.basis.d;
    Coefficients<Lexicographical,T,Index> ret;

    if (Lambda.size()==0) return ret;

    //Determine compression level
    int J=0;        //compression
    if      (d==2) {   J = -std::ceil(2*std::log(tol/((3*tol+3*delta)*kappa))); }
    else if (d==3) {   J = -std::ceil((1./1.5)*std::log(tol/((3*tol+3*delta)*kappa))); }
    else if (d==4) {   J = -std::ceil((1./2.5)*std::log(tol/((3*tol+3*delta)*kappa))); }
    else              { assert(0); }
    std::cerr << "   Estimated compression level for delta=" << delta << " and target tol=" << tol
              << " : " << J << std::endl;

    //Assemble sparse matrix B
    unsigned long N = Lambda.size();
    //std::cerr << "    Assembling of B started with N=" << N << std::endl;

    flens::SparseGeMatrix<CRS<T,CRS_General> > B(N,N);
    A.toFlensSparseMatrix(Lambda,Lambda,B,J);

/*
    A.extendFlensSparseMatrix(supp(w),Extension,sparseMat_A);
    SparseMatrixT B(N,N,sparseMat_A.initializer());
    B.finalize();
*/

    //std::cerr << "    Assembling of B finished." << std::endl;

    //Compute approximate initial residual
    //std::cerr << "    Solving Bx=r0." << std::endl;
    Coefficients<Lexicographical,T,Index> r0, APPLY_Aw;
    APPLY_Aw = A.apply(w, tol/3.);
    r0 = g - P(APPLY_Aw, Lambda);

    DenseVector<Array<T> > rhs(N), x(N), res(N), Bx(N);
    int row_count=1;
    const_coeff_it r0_end = r0.end();
    for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row) {
        const_coeff_it it = r0.find(*row);

        if (it != r0_end) rhs((*row).linearindex) = (*it).second;
        else              rhs((*row).linearindex) = 0.;
    }

    int iters = lawa::cg(B,x,rhs,tol/3.);
    linsolve_iterations.push_back(iters);
    Bx = B*x;
    res= Bx-rhs;
    T lin_res = std::sqrt(res*res);
    std::cerr << "    cg-method needed " << iters << " iterations, res=" << lin_res << std::endl;
    assert(lin_res<tol);


    const_coeff_it w_end = w.end();

    for (const_coeff_it it=w.begin(); it!=w.end(); ++it) {
        ret[(*it).first] = (*it).second + x((*it).first.linearindex);
    }
    for (const_set_it row=Extension.begin(); row!=Extension.end(); ++row) {
        ret[*row] = x((*row).linearindex);
    }


    return ret;
}

}   //namespace lawa
