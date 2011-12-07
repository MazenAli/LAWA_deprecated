namespace lawa {

template <typename TestBasis, typename TrialBasis>
LocalOperator<TestBasis, TrialBasis>::LocalOperator(const TestBasis &_test_basis,
                                                    bool test_withDirichletBC,
                                                    const TrialBasis &_trial_basis,
                                                    bool trial_withDirichletBC, const int _offset)
: test_basis(_test_basis), trial_basis(_trial_basis),
  test_localtransform(test_basis,test_withDirichletBC),
  trial_localtransform(trial_basis,trial_withDirichletBC),
  integral(test_basis,trial_basis), offset(_offset)
{

}

template <typename TestBasis, typename TrialBasis>
void
LocalOperator<TestBasis, TrialBasis>::evalA
(int l, const CoefficientsByLevel<T> &d, const TreeCoefficients1D<T> &c,
        CoefficientsByLevel<T> &PhiPiCheck_vs_v, TreeCoefficients1D<T> &PsiLambdaCheck_vs_v)
{

    if (PhiPiCheck_vs_v.size()==0 && PsiLambdaCheck_vs_v[l].size()==0) return;

    Timer time;
//    time.start();
    CoefficientsByLevel<T> PhiPiCheck2_vs_v;
    if (c[l].size()!=0) {
        const_by_level_it p_cl_end = c[l].end();
        for (const_by_level_it mu=PhiPiCheck_vs_v.begin(); mu!=PhiPiCheck_vs_v.end(); ++mu) {
            long k_col = (*mu).first;
            long k_row_first = std::max(k_col-offset, (long)trial_basis.rangeJ(l).firstIndex());
            long k_row_last  = std::min(k_col+offset, (long)trial_basis.rangeJ(l).lastIndex());
            for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {
                if (overlap(test_basis.mra.phi.support(l,k_col), trial_basis.psi.support(l,k_row)) > 0 ) {
                    const_by_level_it p_entry_cl = c[l].find(k_row);
                    if (p_entry_cl!=p_cl_end) {
                        PhiPiCheck2_vs_v[(*mu).first] = 0.;
                        PhiPiCheck_vs_v.erase((*mu).first);
                        break;
                    }
                }

            }
        }
    }
//    time.stop();
//    T time_PiCheck = time.elapsed();
//    std::cout << "l = " << l << ", PhiPiCheck1_vs_v = " << PhiPiCheck_vs_v << std::endl;
//    std::cout << "l = " << l << ", PhiPiCheck2_vs_v = " << PhiPiCheck2_vs_v << std::endl;

    // Splitting of d
    CoefficientsByLevel<T> d1, d2;
    d1 = d;
//    time.start();
    for (const_by_level_it mu=PsiLambdaCheck_vs_v[l].begin(); mu!=PsiLambdaCheck_vs_v[l].end(); ++mu) {
        long k_col = (*mu).first;
        long k_row_first = std::max(k_col-offset, (long)trial_basis.mra.rangeI(l).firstIndex());
        long k_row_last  = std::min(k_col+offset, (long)trial_basis.mra.rangeI(l).lastIndex());

        for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {
            if (overlap(test_basis.psi.support(l,k_col),
                        trial_basis.mra.phi.support(l,k_row) ) >0 ) {
                const_by_level_it p_entry_d1 = d1.find(k_row);
                if (p_entry_d1!=d1.end()) {
                    d2[k_row] = (*p_entry_d1).second;
                    d1.erase(p_entry_d1);
                }
            }
        }
    }
//    time.stop();
//    T time_d2 = time.elapsed();
//    std::cout << "l = " << l << ", d1 = " << d1 << std::endl;
//    std::cout << "l = " << l << ", d2 = " << d2 << std::endl;

    // Compute underlinePiCheck
//    time.start();
    CoefficientsByLevel<T> PhiPiunderlineCheck_vs_v;
    trial_localtransform.reconstruct(PhiPiCheck2_vs_v, PsiLambdaCheck_vs_v[l], l, PhiPiunderlineCheck_vs_v);
//    time.stop();
//    T time_Piunderline = time.elapsed();
//    std::cout << "l = " << l << ", PhiPiunderlineCheck_vs_v = " << PhiPiunderlineCheck_vs_v  << std::endl;

    //Compute <\Phi_{PiCheck1} , v >
//    time.start();
    const_by_level_it p_d_end = d.end();
    for (by_level_it mu=PhiPiCheck_vs_v.begin(); mu!=PhiPiCheck_vs_v.end(); ++mu) {
        T val = 0.;
        long k_col = (*mu).first;
        long k_row_first = std::max(k_col-offset, (long)trial_basis.mra.rangeI(l).firstIndex());
        long k_row_last  = std::min(k_col+offset, (long)trial_basis.mra.rangeI(l).lastIndex());
        for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {
            const_by_level_it p_entry_d = d.find(k_row);
            if (p_entry_d!=p_d_end) {
                val += integral(l, k_col, XBSpline, 0,
                                l, k_row, XBSpline, 0) * (*p_entry_d).second;
                //val += 0.5 * (*p_entry_d).second;
            }
        }
        (*mu).second = val;
    }
//    time.stop();
//    T time_integration1 = time.elapsed();

//    std::cout << "l = " << l << ", PhiPiCheck1_vs_v = " << PhiPiCheck_vs_v << std::endl;

    // Compute underline d
//    time.start();
    CoefficientsByLevel<T> underline_d;
    trial_localtransform.reconstruct(d2, c[l], l, underline_d);
//    time.stop();
//    T time_underline_d = time.elapsed();
//    std::cout << "l = " << l << ", underline_d = " << underline_d << std::endl;


    //Recursive call of eval
//    time.start();
    this->evalA(l+1, underline_d, c, PhiPiunderlineCheck_vs_v, PsiLambdaCheck_vs_v);
//    time.stop();
//    T time_evalA = time.elapsed();

//    time.start();
    test_localtransform.decompose_(PhiPiunderlineCheck_vs_v, l+1,
                                   PhiPiCheck2_vs_v, PsiLambdaCheck_vs_v[l]);
//    time.stop();
//    T time_decompose = time.elapsed();
//    std::cout << "l = " << l << ", PhiPiCheck2_vs_v = " << PhiPiCheck2_vs_v << std::endl;
//    std::cout << "l = " << l << ", PsiLambdaCheck_vs_v[l] = " << PsiLambdaCheck_vs_v[l] << std::endl;

//    time.start();
    if (PhiPiCheck2_vs_v.size()!=0 && d1.size()!=0) {
        const_by_level_it p_d1_end = d1.end();
        for (by_level_it mu=PhiPiCheck2_vs_v.begin(); mu!=PhiPiCheck2_vs_v.end(); ++mu) {
            T val = 0.;
            long k_col = (*mu).first;
            long k_row_first = std::max(k_col-offset, (long)trial_basis.mra.rangeI(l).firstIndex());
            long k_row_last  = std::min(k_col+offset, (long)trial_basis.mra.rangeI(l).lastIndex());
            for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {
                const_by_level_it p_entry_d1 = d1.find(k_row);
                if (p_entry_d1!=p_d1_end) {
                    val += integral(l, k_col, XBSpline, 0,
                                    l, k_row, XBSpline, 0) * (*p_entry_d1).second;
                    //val += 0.5 * (*p_entry_d1).second;
                }
            }
            (*mu).second += val;
        }
    }
//    time.stop();
//    T time_integration2 = time.elapsed();

    PhiPiCheck_vs_v+= PhiPiCheck2_vs_v;
/*
    std::cerr << "***  l = " << l << "  ***" << std::endl;
    std::cerr << "time_PiCheck =      " << time_PiCheck << std::endl;
    std::cerr << "time_d2 =           " << time_d2 << std::endl;
    std::cerr << "time_Piunderline =  " << time_Piunderline << std::endl;
    std::cerr << "time_integration1 = " << time_integration1 << std::endl;
    std::cerr << "time_underline_d  = " << time_underline_d << std::endl;
//    std::cerr << "time_evalA =        " << time_evalA << std::endl;
    std::cerr << "time_decompose =    " << time_decompose << std::endl;
    std::cerr << "time_integration2 = " << time_integration2 << std::endl;
*/
    return;
}


}   // namespace lawa
