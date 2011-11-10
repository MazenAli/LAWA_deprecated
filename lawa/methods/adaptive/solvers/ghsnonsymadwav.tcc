namespace lawa {

template <typename T, typename Index, typename AdaptiveOperator, typename RHS,
          typename PP_AdaptiveOperator, typename PP_RHS>
GHS_NONSYM_ADWAV<T,Index,AdaptiveOperator,RHS,PP_AdaptiveOperator,PP_RHS>
::GHS_NONSYM_ADWAV(AdaptiveOperator &_A, RHS &_F, PP_AdaptiveOperator &_PP_A, PP_RHS &_PP_F,
                   bool _optimized_grow, bool _assemble_matrix)
    : A(_A), F(_F), PP_A(_PP_A), PP_F(_PP_F),
      optimized_grow(_optimized_grow), assemble_matrix(_assemble_matrix),
      cA(A.cA), CA(A.CA), kappa(A.kappa),
      alpha(0.), omega(0.), gamma(0.), theta(0.), eps(0.)
{
    omega = 0.01;
    alpha = 1./std::sqrt(kappa)-(1.+1./std::sqrt(kappa))*omega-0.00001;
    gamma = 0.5 * (1./6.) * 1./sqrt(kappa) * (alpha-omega)/(1+omega);
    theta = 2./7.;

}

template <typename T, typename Index, typename AdaptiveOperator, typename RHS,
          typename PP_AdaptiveOperator, typename PP_RHS>
Coefficients<Lexicographical,T,Index>
GHS_NONSYM_ADWAV<T,Index,AdaptiveOperator,RHS,PP_AdaptiveOperator,PP_RHS>
::SOLVE(T nuM1, T _eps, const char *filename, int NumOfIterations, T H1norm)
{
    T eps = _eps;
    Coefficients<Lexicographical,T,Index> w_k, w_kP1, g_kP1;
    IndexSet<Index> Lambda_kP1;
    T nu_kM1 = nuM1;
    T nu_k   = 0.;
    T total_time=0.;

    std::cerr << "GHS-NONSYM-ADWAV-SOLVE has started with the following parameters: " << std::endl;
    std::cerr << "  alpha=" << alpha << ", omega=" << omega << ", gamma=" << gamma << ", theta="
              << theta << std::endl;
    std::cerr << "  cA=" << cA << ", CA=" << CA << ", kappa=" << kappa << std::endl;

    std::ofstream file(filename);

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


        Coefficients<Lexicographical,T,Index> PP_Au, PP_f, u;
        u = w_k;
        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            u[(*it).first] *= A.prec((*it).first);
            u[(*it).first] *= 1./PP_A.prec((*it).first);
        }
        PP_f = F(supp(u));
        PP_f = PP_F(supp(u));
        T fu = u*PP_f;
        PP_Au = PP_A.mv(supp(u), u);
        T uAu = u*PP_Au;
        T Error_H_energy = sqrt(fabs(std::pow(H1norm,2.)- 2*fu + uAu));

        Coefficients<Lexicographical,T,Index> Au_M_f;
        A.apply(w_k,0.,Au_M_f,NoTrans);
        Au_M_f -= F(1e-7);
        file << w_k.size() << " " << total_time << " " <<  nu_k << " "
                         << Error_H_energy << " " << Au_M_f.norm(2.) << std::endl;

        /*
        std::stringstream coeff_filename;
        coeff_filename << "adwav_coeff_" << w_k.size();
        Coefficients<AbsoluteValue,T,Index1D> u_abs;
        u_abs = w_k;
        plotCoeff(u_abs, A.basis, coeff_filename.str().c_str());
        */
        /*
        std::stringstream coefffile;
        coefffile << "s_adwav_coeffs_" << i;
        plotScatterCoeff2D(w_k, A.basis.first, A.basis.second, coefffile.str().c_str());
        */

        time.start();
        std::cerr << "   GALSOLVE started with #Lambda = " << Lambda_kP1.size()  << std::endl;

        /*
         * Attention: update in rhs1d is set! Linear complexity still holds with better convergence!
         */
        Coefficients<Lexicographical,T,Index> rhs;
        //g_kP1 = F(Lambda_kP1);                // update rhs vector
        rhs = F(gamma*nu_k);
        g_kP1 = P(rhs,Lambda_kP1);  // compute with restriction, otherwise GROW does not work
        w_kP1 = this->GALSOLVE(Lambda_kP1, Extension, g_kP1, w_k, (1+gamma)*nu_k, gamma*nu_k);

        std::cerr << "  GALSOLVE finished." << std::endl;
/*
        std::stringstream coefffile;
        coefffile << "ghs_adwav_coeffs_" << i;
        plotScatterCoeff2D(w_kP1, A.basis.first, A.basis.second, coefffile.str().c_str());
*/
        nu_kM1 = nu_k;
        w_k = w_kP1;
        time.stop();
        total_time += time.elapsed();
        std::cerr << std::endl;
    }
    return w_k;
}

template <typename T, typename Index, typename AdaptiveOperator, typename RHS,
          typename PP_AdaptiveOperator, typename PP_RHS>
IndexSet<Index>
GHS_NONSYM_ADWAV<T,Index,AdaptiveOperator,RHS,PP_AdaptiveOperator,PP_RHS>
::GROW(const Coefficients<Lexicographical,T,Index> &w, T nu_bar, T &nu)
{
    T zeta = 2.*(omega*nu_bar)/(1-omega);
    T r_norm = 0.;
    Coefficients<Lexicographical,T,Index> r, rhs, help;
    while (1) {
        zeta *= 0.5;

        std::cerr << "GROW: zeta = " << zeta << std::endl;
        A.apply(w, 0.25*zeta, help, cxxblas::NoTrans);
        help = THRESH(help,0.1*zeta);
        A.apply(help, 0.25*zeta, r, cxxblas::Trans);
        std::cerr << "      size of r: " << r.size() << std::endl;
        help = F(0.25*zeta);
        A.apply(help, 0.25*zeta, rhs, cxxblas::Trans);
        r -= rhs;

        r_norm = r.norm(2.);
        nu = r_norm + zeta;
        std::cerr << "    zeta = " << zeta << ", r_norm = " << r_norm << ", omega*r_norm = "
                  << omega*r_norm << ", nu = " << nu << std::endl;
        if (nu <= eps) break;
        if (zeta<=omega*r_norm) break;
    }

    IndexSet<Index> Lambda;
    long double P_Lambda_r_norm_square = 0.0L;
    if (w.size() > 0) {
        for (const_coeff_it it=w.begin(); it!=w.end(); ++it) {
            //Lambda.insert((*it).first); -> Extension is calculated!!
            P_Lambda_r_norm_square += std::pow(r[(*it).first],2.);
            r.erase((*it).first);
        }
    }

    T threshbound = std::sqrt(1-alpha*alpha) * r.norm(2.)/std::sqrt(T(r.size()));
    Coefficients<Bucket,T,Index> r_bucket;
    /*
    Coefficients<AbsoluteValue,T,Index> r_abs;
    r_abs = r;
    int count=1;
    for (const_coeff_abs_it it=r_abs.begin(); it!=r_abs.end(); ++it,++count) {
        std::cerr << "   " << (*it).first << ", " << (*it).second << std::endl;
        if (count >10) break;
    }
    getchar();
    */

    r_bucket.bucketsort(r, threshbound);
    //std::cerr << r_bucket << std::endl;

    std::cout << "  GROW: size of r = " << r.size() << std::endl;
    std::cerr << "   Threshbound = " << threshbound << std::endl;
    std::cerr << "   Before extension: ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_r_norm_square)
              << ", alpha*r_norm = " << alpha*r_norm << std::endl;
    if (nu > eps) {
        int currentDoF = w.size();
        for (int i=0; i<(int)r_bucket.bucket_ell2norms.size(); ++i) {
            P_Lambda_r_norm_square += std::pow(r_bucket.bucket_ell2norms[i],2.0L);
            std::cerr << "   ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_r_norm_square) << std::endl;
            int addDoF = r_bucket.addBucketToIndexSet(Lambda,i,currentDoF);
            currentDoF += addDoF;
            if (P_Lambda_r_norm_square >= alpha*r_norm*alpha*r_norm) {
                if (optimized_grow) {
                    int addDoF = r_bucket.addBucketToIndexSet(Lambda,i+1,currentDoF);
                    currentDoF += addDoF;
                }
                break;
            }
        }
    }
    //getchar();
    return Lambda;
}

template <typename T, typename Index, typename AdaptiveOperator, typename RHS,
          typename PP_AdaptiveOperator, typename PP_RHS>
Coefficients<Lexicographical,T,Index>
GHS_NONSYM_ADWAV<T,Index,AdaptiveOperator,RHS,PP_AdaptiveOperator,PP_RHS>
::GALSOLVE(const IndexSet<Index> &Lambda, const IndexSet<Index> &Extension,
           const Coefficients<Lexicographical,T,Index> &g, const Coefficients<Lexicographical,T,Index> &w,
           T delta, T tol)
{

    if (assemble_matrix) {
        IndexSet<Index> LambdaCol, LambdaRow;
        LambdaCol = Lambda + Extension;
        Coefficients<Lexicographical,T,Index> help1, help2, help3;
        for (const_set_it it=LambdaCol.begin(); it!=LambdaCol.end(); ++it) {
            help1[*it] = 1.;
        }
        A.apply(help1,0.,help2,NoTrans);
        LambdaRow = supp(help2);

    //    T res1 = 0.;
    //    int n=CGLS_Solve(LambdaRow, LambdaCol, A, help3, g, res1, tol/3., 1000);
    //    return help3;

        flens::SparseGeMatrix<CRS<T,CRS_General> > B(LambdaRow.size(),LambdaCol.size());
        A.toFlensSparseMatrix(LambdaRow,LambdaCol,B,1,false);

        std::cerr << "LambdaRow.size() = " << LambdaRow.size() << ", LambdaCol.size() = " << LambdaCol.size() << std::endl;

        DenseVector<Array<T> > g_vec(LambdaRow.size()), x_vec(LambdaCol.size()), Ax_vec(LambdaRow.size());
        const_coeff_it g_end = g.end();
        const_coeff_it w_end = w.end();
        int row_count=1;
        for (const_set_it row=LambdaRow.begin(); row!=LambdaRow.end(); ++row) {
            const_coeff_it it = g.find(*row);
            if (it != g_end) g_vec(row_count) = (*it).second;
            else             g_vec(row_count) = 0.;
            ++row_count;
        }
        row_count=1;
        for (const_set_it row=LambdaCol.begin(); row!=LambdaCol.end(); ++row) {
            const_coeff_it it2 = w.find(*row);
            if (it2 != w_end) x_vec(row_count) = (*it2).second;
            else              x_vec(row_count) = 0.;
            ++row_count;
        }
        int iters = lawa::cgls(B,x_vec,g_vec,tol);

        linsolve_iterations.push_back(iters);

        Coefficients<Lexicographical,T,Index> x;
        row_count=1;
        for (const_set_it row=LambdaCol.begin(); row!=LambdaCol.end(); ++row) {
            x[(*row)] = x_vec(row_count);
            ++row_count;
        }
        return x;
    }
    else {
        std::cerr << "     Not yet implemented. Sorry..." << std::endl;
        exit(1);
    }



    /*


    T alpha, beta, gammaPrev, gamma, b_norm;
    Coefficients<Lexicographical,T,Index1D> b,x2;
    Coefficients<Lexicographical,T,Index1D> r2, q, s, p;
    b = g;
    x2 = w;
    b_norm = b.norm(2.);
    A.apply(x2,0.,r2,NoTrans);
    r2 -= b;
    r2 *= -1.;
    std::cerr << "  r = " << r2 << std::endl;
    A.apply(r2,0.,LambdaCol,s,Trans);
    p = s;
    gammaPrev = s*s;

    for (int k=1; k<=50; k++) {
       A.apply(p,0.,q,NoTrans);   //q = A*p;
       alpha = gammaPrev/(q*q);
       x2 +=   alpha *p;
       r2 -=   alpha*q;
       A.apply(r2,0.,LambdaCol,s,Trans);  // flens::blas::mv(cxxblas::Trans, typename _cg<VB>::T(1), A, r, typename _cg<VB>::T(0), s);

       gamma = s*s;
       //std::cout << "Iteration " << k << ": gamma = " << gamma << std::endl;

       if (sqrt(gamma)<=tol) {
           std::cerr << "   inner iterations: " << k << std::endl;
           break;
       }

       beta  = gamma/gammaPrev;
       p *= beta;
       p += s;
       gammaPrev = gamma;
    }

    std::cerr << "   x2 = " << x2 << std::endl;


    Coefficients<Lexicographical,T,Index> Ax, AtAx, Atg, residual;
    A.apply(x,0.,Ax,NoTrans);
    A.apply(Ax,0.,LambdaCol,AtAx,Trans);
    A.apply(g,0.,LambdaCol,Atg,Trans);
    residual = AtAx - Atg;
    std::cerr << "   ||At(Ax-g)||_2 = " << residual.norm(2.) << std::endl;

    Coefficients<Lexicographical,T,Index> Ax2, AtAx2, residual2;
    A.apply(x2,0.,Ax2,NoTrans);
    A.apply(Ax2,0.,LambdaCol,AtAx2,Trans);
    residual2 = AtAx2 - Atg;
    std::cerr << "   ||At(Ax2-g)||_2 = " << residual2.norm(2.) << std::endl;

    */


}

}   // namespace lawa
