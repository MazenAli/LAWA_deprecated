namespace lawa {

template <typename T, typename PrimalTestBasis, typename PrimalTrialBasis>
void
eval(int l, const PrimalTestBasis &test_basis, const PrimalTrialBasis &trial_basis,
            std::map<int, Coefficients<Lexicographical,T,Index1D> >          &c,
            const Coefficients<Lexicographical,T,Index1D>                    &d_lM1,
            std::map<int, IndexSet<Index1D> >                                &Lhd,
            const IndexSet<Index1D>                                          &SquareCap_lM1,
            Coefficients<Lexicographical,T,Index1D>                          &upsilon_vs_v_lM1_SquareCap_lM1,
            Coefficients<Lexicographical,T,Index1D>                          &theta_vs_v_lM1_Lhd)
{
    typedef IndexSet<Index1D>::const_iterator                                  const_set1d_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator   const_coeff1d_it;
    typedef Integral<Gauss, PrimalTestBasis, PrimalTrialBasis>                 Integral_Psi_Psi;

    //std::cerr << "l = " << l << std::endl << std::endl;

    if (SquareCap_lM1.size()==0 && Lhd.count(l)==0) return;

    // Splitting of SquareCap_lM1
    IndexSet<Index1D> SquareCap_lM1_1, SquareCap_lM1_2;
    SquareCap_lM1_1 = SquareCap_lM1;
    if (c.count(l)>0) {
        for (const_set1d_it mu=SquareCap_lM1.begin(); mu!=SquareCap_lM1.end(); ++mu) {
            //std::cerr << *mu << " : " << test_basis.mra.phi.support((*mu).j,(*mu).k) << std::endl;
            for (const_coeff1d_it lambda=c[l].begin(); lambda!=c[l].end(); ++lambda) {
                //std::cerr << (*lambda).first << " : " << trial_basis.psi.support((*lambda).first.j,(*lambda).first.k) << std::endl;
                if (overlap(test_basis.mra.phi.support((*mu).j,(*mu).k),
                            trial_basis.psi.support((*lambda).first.j,(*lambda).first.k))>0 ) {
                    SquareCap_lM1_2.insert(*mu);
                    SquareCap_lM1_1.erase(*mu);
                    break;
                }
            }
        }
    }
//    std::cout << "l = " << l << ", PiCheck1 = " << SquareCap_lM1_1 << std::endl;
//    std::cout << "l = " << l << ", PiCheck2 = " << SquareCap_lM1_2 << std::endl;

    // Splitting of d_lM1
    Coefficients<Lexicographical,T,Index1D> d_lM1_1, d_lM1_2;
    d_lM1_1 = d_lM1;
    if (Lhd.count(l)>0) {
        for (const_coeff1d_it mu=d_lM1.begin(); mu!=d_lM1.end(); ++mu) {
            for (const_set1d_it lambda=Lhd[l].begin(); lambda!=Lhd[l].end(); ++lambda) {
                if (overlap(trial_basis.mra.phi.support((*mu).first.j,(*mu).first.k),
                            test_basis.psi.support((*lambda).j,(*lambda).k))>0 ) {
                    d_lM1_2[(*mu).first] = (*mu).second;
                    d_lM1_1.erase((*mu).first);
                    break;
                }
            }
        }
    }
    std::cout << "l = " << l << ", d1 = " << d_lM1_1 << std::endl;
    std::cout << "l = " << l << ", d2 = " << d_lM1_2 << std::endl;

    //Compute SquareCap_l
    IndexSet<Index1D> SquareCap_lM1_2_u_Lhd_l, SquareCap_l;
    SquareCap_lM1_2_u_Lhd_l = SquareCap_lM1_2;
    if (Lhd.count(l)>0) {   SquareCap_lM1_2_u_Lhd_l += Lhd[l];  }
    if (SquareCap_lM1_2_u_Lhd_l.size()>0) {
        computeLocalReconstruction<T,PrimalTestBasis>(SquareCap_lM1_2_u_Lhd_l, test_basis, l, SquareCap_l);
    }
    //std::cout << SquareCap_lM1_2_u_Lhd_l  << std::endl;
//    std::cout << "l = " << l << ", underlinePiCheck = " << SquareCap_l  << std::endl;


    //Compute <\upsilon_{\mu} , v_{l-1} > for \mu \in \SquareCap_{l-1}^1
    Coefficients<Lexicographical,T,Index1D> upsilon_vs_v_lM1_SquareCap_lM1_1;
    Integral_Psi_Psi integral(test_basis,trial_basis);
    for (const_set1d_it mu=SquareCap_lM1_1.begin(); mu!=SquareCap_lM1_1.end(); ++mu) {
        T val = 0.;
        for (const_coeff1d_it lambda=d_lM1.begin(); lambda!=d_lM1.end(); ++lambda) {
            val += integral((*mu).j,(*mu).k,(*mu).xtype,0,
                             (*lambda).first.j,(*lambda).first.k,(*lambda).first.xtype,0) * (*lambda).second;
        }
        upsilon_vs_v_lM1_SquareCap_lM1_1[*mu] = val;
    }
//    std::cout << "l = " << l << ", PhiPiCheck1_vs_v = "  << upsilon_vs_v_lM1_SquareCap_lM1_1 << std::endl;


    //Compute lsf e
    Coefficients<Lexicographical,T,Index1D> e_multi, e, tmp;
    e_multi = d_lM1_2;
    if (c.count(l)>0)   { e_multi+= c[l]; }
    Timer time;
    time.start();
    computeLocalReconstruction(e_multi, trial_basis, l, e);
    time.stop();
    T time_oldLocalReconstruction = time.elapsed();
    LocalRefinement<PrimalTrialBasis> LocalRefine(trial_basis,false);
    time.start();
    LocalRefine.reconstruct(e_multi, l, tmp);
    time.stop();
    tmp -= e;
    //std::cerr << "Testing local refinement: error = " << tmp.norm(2.)
    //          << ", time was " << time.elapsed() << "(" << time_oldLocalReconstruction << ")"<< std::endl;

    //std::cout << "e_multi = " << e_multi << std::endl;
//    std::cout << "l = " << l << ", underline_d = " << e << std::endl;

    //Recursive call of eval
    Coefficients<Lexicographical,T,Index1D> UpsilonVsV_l;
    eval<T, PrimalTestBasis, PrimalTrialBasis>(l+1, test_basis, trial_basis, c, e, Lhd, SquareCap_l, UpsilonVsV_l, theta_vs_v_lM1_Lhd);

    Coefficients<Lexicographical,T,Index1D> upsilon_vs_v_lM1_SquareCap_lM1_2;
    //Coefficients<Lexicographical,T,Index1D> theta_vs_v_lM1_Lhd_l;
    //std::cerr << "l = " << l << ", current Lambda = " << SquareCap_lM1_2+Lhd[l] << std::endl;
    Coefficients<Lexicographical,T,Index1D> test_upsilon_vs_v_lM1_SquareCap_lM1_2, test_theta_vs_v_lM1_Lhd;
    test_theta_vs_v_lM1_Lhd = theta_vs_v_lM1_Lhd;
//    std::cerr << "l = " << l << ", PhiPiunderlineCheck_vs_v = " << UpsilonVsV_l << std::endl;
    time.start();
    computeLocalDecomposition_(UpsilonVsV_l, test_basis, l+1, SquareCap_lM1_2+Lhd[l],
                               upsilon_vs_v_lM1_SquareCap_lM1_2, theta_vs_v_lM1_Lhd);
    time.stop();
//    std::cout << "l = " << l << ", PhiPiCheck2_vs_v = " << upsilon_vs_v_lM1_SquareCap_lM1_2 << std::endl;
//    std::cout << "l = " << l << ", PsiLambdaCheck_vs_v[l] = " << theta_vs_v_lM1_Lhd << std::endl;

    T time_oldLocalDecompose = time.elapsed();
    LocalRefinement<PrimalTestBasis> LocalRefineTest(test_basis,false);
    time.start();
    LocalRefineTest.decompose_(UpsilonVsV_l, SquareCap_lM1_2+Lhd[l], test_upsilon_vs_v_lM1_SquareCap_lM1_2, test_theta_vs_v_lM1_Lhd);
    time.stop();
    test_upsilon_vs_v_lM1_SquareCap_lM1_2 -= upsilon_vs_v_lM1_SquareCap_lM1_2;
    test_theta_vs_v_lM1_Lhd -= theta_vs_v_lM1_Lhd;
    //std::cerr << "Testing local decompose: error = " << test_upsilon_vs_v_lM1_SquareCap_lM1_2.norm(2.)
    //          << " , " << test_theta_vs_v_lM1_Lhd.norm(2.)
    //          << ", time was " << time.elapsed() << "(" << time_oldLocalDecompose << ")"<< std::endl;

    //std::cerr << "l = " << l << ", theta_vs_v_lM1_Lhd " << theta_vs_v_lM1_Lhd << std::endl;
    //std::cerr << "l = " << l << ", upsilon_vs_v_lM1_SquareCap_lM1_2 " << upsilon_vs_v_lM1_SquareCap_lM1_2 << std::endl;
    for (const_set1d_it mu=SquareCap_lM1_2.begin(); mu!=SquareCap_lM1_2.end(); ++mu) {
        T val = 0.;
        for (const_coeff1d_it lambda=d_lM1_1.begin(); lambda!=d_lM1_1.end(); ++lambda) {
             val += integral((*mu).j,(*mu).k,(*mu).xtype,0,
                             (*lambda).first.j,(*lambda).first.k,(*lambda).first.xtype,0) * (*lambda).second;
        }
        upsilon_vs_v_lM1_SquareCap_lM1_2[*mu] += val;
    }
    //std::cout << "l = " << l << ", upsilon_vs_v_lM1_SquareCap_lM1_1 = " << upsilon_vs_v_lM1_SquareCap_lM1_1 << std::endl;
    //std::cerr << "l = " << l << ", upsilon_vs_v_lM1_SquareCap_lM1_2 = " << upsilon_vs_v_lM1_SquareCap_lM1_2 << std::endl;

    upsilon_vs_v_lM1_SquareCap_lM1 = upsilon_vs_v_lM1_SquareCap_lM1_1;
    upsilon_vs_v_lM1_SquareCap_lM1+= upsilon_vs_v_lM1_SquareCap_lM1_2;
    //theta_vs_v_lM1_Lhd  += theta_vs_v_lM1_Lhd_l;
    //std::cerr << "l = " << l << ", upsilon_vs_v_lM1_SquareCap_lM1 " << upsilon_vs_v_lM1_SquareCap_lM1 << std::endl;
    //std::cerr << "l = " << l << ", theta_vs_v_lM1_Lhd_l " << theta_vs_v_lM1_Lhd_l << std::endl;
    return;

}

}   // namespace lawa
