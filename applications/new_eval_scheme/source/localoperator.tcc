namespace lawa {

template <typename TestBasis, typename TrialBasis, typename BilinearForm, typename Preconditioner>
LocalOperator<TestBasis, TrialBasis, BilinearForm, Preconditioner>
::LocalOperator(const TestBasis &_test_basis, bool test_withDirichletBC,
                const TrialBasis &_trial_basis, bool trial_withDirichletBC, const int _offset,
                BilinearForm &_Bil, const Preconditioner &_Prec)
: test_basis(_test_basis), trial_basis(_trial_basis),
  Bil(_Bil), Prec(_Prec),
  test_localtransform(test_basis,test_withDirichletBC),
  trial_localtransform(trial_basis,trial_withDirichletBC),
  offset(_offset)
{

}

template <typename TestBasis, typename TrialBasis, typename BilinearForm, typename Preconditioner>
void
LocalOperator<TestBasis, TrialBasis, BilinearForm, Preconditioner>
::scale_wrt_trialbasis(const TreeCoefficients1D<T> &x, TreeCoefficients1D<T> &y)
{
    for (const_by_level_it it=x.bylevel[trial_basis.j0-1].map.begin(); it!=x.bylevel[trial_basis.j0-1].map.end(); ++it) {
        y[trial_basis.j0-1].map[((*it).first)] = (*it).second * Prec(XBSpline,trial_basis.j0,(*it).first);
    }
    for (int l=trial_basis.j0; l<=JMAX; ++l) {
        if (x.bylevel[l].map.size()!=0) {
            for (const_by_level_it it=x.bylevel[l].map.begin(); it!=x.bylevel[l].map.end(); ++it) {
                y[l].map[((*it).first)] = (*it).second * Prec(XWavelet,l,(*it).first);
            }
        }
        else return;
    }
}

template <typename TestBasis, typename TrialBasis, typename BilinearForm, typename Preconditioner>
void
LocalOperator<TestBasis, TrialBasis, BilinearForm, Preconditioner>
::computePhiPi2(int l, const CoefficientsByLevel<T> &cl, CoefficientsByLevel<T> &PhiPiCheck_vs_v,
                CoefficientsByLevel<T> &PhiPiCheck2_vs_v) /*const*/
{
    if (cl.map.size()==0) return;
    if (cl.map.size()>PhiPiCheck_vs_v.map.size()) {
        const_by_level_it p_cl_end = cl.map.end();
        for (const_by_level_it mu=PhiPiCheck_vs_v.map.begin(); mu!=PhiPiCheck_vs_v.map.end(); ++mu) {
            long k_col = (*mu).first;
            long k_row_first = std::max(k_col-offset, (long)trial_basis.rangeJ(l).firstIndex());
            long k_row_last  = std::min(k_col+offset, (long)trial_basis.rangeJ(l).lastIndex());
            for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {
                if (overlap(test_basis.mra.phi.support(l,k_col), trial_basis.psi.support(l,k_row)) > 0 ) {
                    const_by_level_it p_entry_cl = cl.map.find(k_row);
                    if (p_entry_cl!=p_cl_end) {
                        PhiPiCheck2_vs_v.map[(*mu).first] = 0.;
                        PhiPiCheck_vs_v.map.erase((*mu).first);
                        break;
                    }
                }
            }
        }
    }
    else {
        for (const_by_level_it mu=cl.map.begin(); mu!=cl.map.end(); ++mu) {
            long k_col = (*mu).first;
            long k_row_first = std::max(k_col-offset, (long)test_basis.mra.rangeI(l).firstIndex());
            long k_row_last  = std::min(k_col+offset, (long)test_basis.mra.rangeI(l).lastIndex());
            for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {
                if (overlap(test_basis.mra.phi.support(l,k_row), trial_basis.psi.support(l,k_col)) > 0 ) {
                    const_by_level_it p_PhiPiCheck_vs_v = PhiPiCheck_vs_v.map.find(k_row);
                    if (p_PhiPiCheck_vs_v!=PhiPiCheck_vs_v.map.end()) {
                        PhiPiCheck2_vs_v.map[k_row] = 0.;
                        PhiPiCheck_vs_v.map.erase(k_row);
                    }
                }
            }
        }
    }
}

template <typename TestBasis, typename TrialBasis, typename BilinearForm, typename Preconditioner>
void
LocalOperator<TestBasis, TrialBasis, BilinearForm, Preconditioner>
::computed2(int l, const CoefficientsByLevel<T> &PsiLambdaCheck_vs_v_l,
            const CoefficientsByLevel<T> &d, CoefficientsByLevel<T> &d1, CoefficientsByLevel<T> &d2)
/*const*/
{
    if (PsiLambdaCheck_vs_v_l.map.size()==0) return;
    if (d.map.size()>PsiLambdaCheck_vs_v_l.map.size()) {
        for (const_by_level_it mu=PsiLambdaCheck_vs_v_l.map.begin(); mu!=PsiLambdaCheck_vs_v_l.map.end(); ++mu) {
            long k_col = (*mu).first;
            long k_row_first = std::max(k_col-offset, (long)trial_basis.mra.rangeI(l).firstIndex());
            long k_row_last  = std::min(k_col+offset, (long)trial_basis.mra.rangeI(l).lastIndex());

            for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {
                if (overlap(test_basis.psi.support(l,k_col),
                            trial_basis.mra.phi.support(l,k_row) ) >0 ) {
                    by_level_it p_entry_d1 = d1.map.find(k_row);
                    if (p_entry_d1!=d1.map.end()) {
                        d2.map[k_row] = (*p_entry_d1).second;
                        d1.map.erase(p_entry_d1);
                    }
                }
            }
        }
    }
    else {
        const_by_level_it p_PsiLambdaCheck_vs_v_l_end = PsiLambdaCheck_vs_v_l.map.end();
        for (const_by_level_it mu=d.map.begin(); mu!=d.map.end(); ++mu) {
            long k_col = (*mu).first;
            long k_row_first = std::max(k_col-offset, (long)test_basis.rangeJ(l).firstIndex());
            long k_row_last  = std::min(k_col+offset, (long)test_basis.rangeJ(l).lastIndex());

            for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {
                if (overlap(test_basis.mra.phi.support(l,k_col), trial_basis.psi.support(l,k_row) ) >0 ) {
                    const_by_level_it p_PsiLambdaCheck_vs_v_l = PsiLambdaCheck_vs_v_l.map.find(k_row);
                    if (p_PsiLambdaCheck_vs_v_l!=p_PsiLambdaCheck_vs_v_l_end) {
                        d2.map[k_col] = (*mu).second;
                        d1.map.erase(k_col);
                        break;
                    }
                }
            }
        }
    }
    /*
     * Here seems to be an error...
    else {
        const_by_level_it p_d_end = d.map.end();
        for (const_by_level_it mu=d.map.begin(); mu!=d.map.end(); ++mu) {
            long k_col = (*mu).first;
            long k_row_first = std::max(k_col-offset, (long)test_basis.mra.rangeI(l).firstIndex());
            long k_row_last  = std::min(k_col+offset, (long)test_basis.mra.rangeI(l).lastIndex());

            for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {
                if (overlap(test_basis.psi.support(l,k_col), trial_basis.mra.phi.support(l,k_row) ) >0 ) {
                    const_by_level_it p_d = d.map.find(k_row);
                    if (p_d!=p_d_end) {
                        d2.map[k_col] = (*mu).second;
                        d1.map.erase(k_col);
                        break;
                    }
                }
            }
        }
    }
    */
}

template <typename TestBasis, typename TrialBasis, typename BilinearForm, typename Preconditioner>
void
LocalOperator<TestBasis, TrialBasis, BilinearForm, Preconditioner>
::applyBilinearForm(int l, const CoefficientsByLevel<T> &d, CoefficientsByLevel<T> &PhiPiCheck_vs_v) /*const*/
{
    if (d.map.size()==0 && PhiPiCheck_vs_v.map.size()==0) return;
    if (d.map.size() > PhiPiCheck_vs_v.map.size()) {
        const_by_level_it p_d_end = d.map.end();
        long leftrangebound =  (long)trial_basis.mra.rangeI(l).firstIndex();
        long rightrangebound = (long)trial_basis.mra.rangeI(l).lastIndex();
        for (by_level_it mu=PhiPiCheck_vs_v.map.begin(); mu!=PhiPiCheck_vs_v.map.end(); ++mu) {
            long double val = 0.;
            long k_col = (*mu).first;
            long k_row_first = std::max(k_col-offset, leftrangebound);
            long k_row_last  = std::min(k_col+offset, rightrangebound);
            for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {
                const_by_level_it p_entry_d = d.map.find(k_row);
                if (p_entry_d!=p_d_end) {
                    val += (long double)(Bil(XBSpline, l, k_col, XBSpline, l, k_row) * (*p_entry_d).second);
                }
            }
            (*mu).second += val;
        }
    }

    else {
        const_by_level_it p_PhiPiCheck_vs_v_end = PhiPiCheck_vs_v.map.end();
        long leftrangebound =  (long)test_basis.mra.rangeI(l).firstIndex();
        long rightrangebound = (long)test_basis.mra.rangeI(l).lastIndex();
        for (const_by_level_it mu=d.map.begin(); mu!=d.map.end(); ++mu) {
            long double val = 0.;
            long k_col = (*mu).first;
            long k_row_first = std::max(k_col-offset, leftrangebound);
            long k_row_last  = std::min(k_col+offset, rightrangebound);
            for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {
                by_level_it p_PhiPiCheck_vs_v=PhiPiCheck_vs_v.map.find(k_row);
                if (p_PhiPiCheck_vs_v!=p_PhiPiCheck_vs_v_end) {
                    (*p_PhiPiCheck_vs_v).second
                    += (long double)(Bil(XBSpline, l, k_col, XBSpline, l, k_row) * (*mu).second);
                }
            }
        }
    }

}

template <typename TestBasis, typename TrialBasis, typename BilinearForm, typename Preconditioner>
void
LocalOperator<TestBasis, TrialBasis, BilinearForm, Preconditioner>
::evalA(int l, const CoefficientsByLevel<T> &d, const TreeCoefficients1D<T> &c,
        CoefficientsByLevel<T> &PhiPiCheck_vs_v, TreeCoefficients1D<T> &PsiLambdaCheck_vs_v,
        bool pre_apply_prec) /*const*/
{
    //if (PhiPiCheck_vs_v.map.size()==0 && PsiLambdaCheck_vs_v[l].map.size()==0) return;
    if (d.map.size()==0 && c[l].map.size()==0) return;

    size_t hm_size = COEFFBYLEVELSIZE;
    Timer time;
//    time.start();
    CoefficientsByLevel<T> PhiPiCheck2_vs_v(l,hm_size);
    computePhiPi2(l, c[l], PhiPiCheck_vs_v, PhiPiCheck2_vs_v);
//    time.stop();
//    T time_PiCheck = time.elapsed();

    // Splitting of d
    CoefficientsByLevel<T> d1(l,hm_size), d2(l,hm_size);
    d1 = d;
//    time.start();
    computed2(l, PsiLambdaCheck_vs_v[l], d, d1, d2);
//    time.stop();
//    T time_d2 = time.elapsed();

    // Compute underlinePiCheck
//    time.start();
    CoefficientsByLevel<T> PhiPiunderlineCheck_vs_v(l,hm_size);
    trial_localtransform.reconstruct(PhiPiCheck2_vs_v, PsiLambdaCheck_vs_v[l], l, PhiPiunderlineCheck_vs_v);
//    time.stop();
//    T time_Piunderline = time.elapsed();

    //Compute <\Phi_{PiCheck1} , v >
//    time.start();
    applyBilinearForm(l, d, PhiPiCheck_vs_v);
//    time.stop();
//    T time_integration1 = time.elapsed();

    // Compute underline d
//    time.start();
    CoefficientsByLevel<T> underline_d(l+1,hm_size);
    trial_localtransform.reconstruct(d2, c[l], l, underline_d);
//    time.stop();
//    T time_underline_d = time.elapsed();

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

//    time.start();
    if (PhiPiCheck2_vs_v.map.size()!=0 && d1.map.size()!=0) {
        applyBilinearForm(l, d1, PhiPiCheck2_vs_v);
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

template <typename TestBasis, typename TrialBasis, typename BilinearForm, typename Preconditioner>
void
LocalOperator<TestBasis, TrialBasis, BilinearForm, Preconditioner>
::evalU(int l, const CoefficientsByLevel<T> &d, const TreeCoefficients1D<T> &c,
        CoefficientsByLevel<T> &PhiPiCheck_vs_v, TreeCoefficients1D<T> &PsiLambdaCheck_vs_v,
        bool pre_apply_prec) /*const*/
{
    if (PhiPiCheck_vs_v.map.size()==0 && PsiLambdaCheck_vs_v[l].map.size()==0) return;
    if (d.map.size()==0 && c[l].map.size()==0) return;

//    if (l==trial_basis.j0 && c.maxTreeLevel<ref_j) {
        //lowCostevalU(d, c, PhiPiCheck_vs_v, PsiLambdaCheck_vs_v);
//    }

    size_t hm_size_PhiPiCheck_vs_v = COEFFBYLEVELSIZE;//PhiPiCheck_vs_v.map.size()*4;
    size_t hm_size_d               = COEFFBYLEVELSIZE;//d.map.size()*4;
//    Timer time;
//    time.start();
    CoefficientsByLevel<T> PhiPiCheck2_vs_v(l,hm_size_PhiPiCheck_vs_v);
    computePhiPi2(l, c[l], PhiPiCheck_vs_v, PhiPiCheck2_vs_v);
//    time.stop();
//    T time_PiCheck = time.elapsed();

    //Compute <\Phi_{PiCheck1} , v >
//    time.start();
    applyBilinearForm(l, d, PhiPiCheck_vs_v);
//    time.stop();
//    T time_integration1 = time.elapsed();


    // Compute underlinePiCheck
//    time.start();
    //CoefficientsByLevel<T> PhiPiunderlineCheck_vs_v;
    CoefficientsByLevel<T> PhiPiunderlineCheck_vs_v(l,hm_size_PhiPiCheck_vs_v);
    trial_localtransform.reconstruct(PhiPiCheck2_vs_v, PsiLambdaCheck_vs_v[l], l, PhiPiunderlineCheck_vs_v);
//    time.stop();
//    T time_Piunderline = time.elapsed();

    // Compute underline d
//    time.start();
    CoefficientsByLevel<T> underline_d(l,hm_size_d), d2(l,hm_size_d);
    trial_localtransform.reconstruct(d2, c[l], l, underline_d);
//    time.stop();
//    T time_underline_d = time.elapsed();

    //Recursive call of eval
//    time.start();
    this->evalU(l+1, underline_d, c, PhiPiunderlineCheck_vs_v, PsiLambdaCheck_vs_v);
//    time.stop();
//    T time_evalA = time.elapsed();

//    time.start();
    test_localtransform.decompose_(PhiPiunderlineCheck_vs_v, l+1,
                                   PhiPiCheck2_vs_v, PsiLambdaCheck_vs_v[l]);
//    time.stop();
//    T time_decompose = time.elapsed();

//    time.start();
    if (!(PhiPiCheck2_vs_v.map.size()==0 || d.map.size()==0)) {
        applyBilinearForm(l, d, PhiPiCheck2_vs_v);
    }

//    time.stop();
//    T time_integration2 = time.elapsed();

    PhiPiCheck_vs_v+= PhiPiCheck2_vs_v;

/*
    std::cerr << "***  l = " << l << "  ***" << std::endl;
    std::cerr << "time_PiCheck =      " << time_PiCheck << std::endl;
    std::cerr << "time_Piunderline =  " << time_Piunderline << std::endl;
    std::cerr << "time_integration1 = " << time_integration1 << std::endl;
    std::cerr << "time_underline_d  = " << time_underline_d << std::endl;
//    std::cerr << "time_evalA =        " << time_evalA << std::endl;
    std::cerr << "time_decompose =    " << time_decompose << std::endl;
    std::cerr << "time_integration2 = " << time_integration2 << std::endl;
*/
    return;
}

template <typename TestBasis, typename TrialBasis, typename BilinearForm, typename Preconditioner>
void
LocalOperator<TestBasis, TrialBasis, BilinearForm, Preconditioner>
::evalL(int l, const CoefficientsByLevel<T> &d, const TreeCoefficients1D<T> &c,
        TreeCoefficients1D<T> &PsiLambdaCheck_vs_v, bool pre_apply_prec) /*const*/
{
    if (PsiLambdaCheck_vs_v[l].map.size()==0) return;
    if (d.map.size()==0 && c[l].map.size()==0) return;

    size_t hm_size = COEFFBYLEVELSIZE;
    Timer time;

    // Splitting of d
    CoefficientsByLevel<T> d1(l,hm_size), d2(l,hm_size);
    //CoefficientsByLevel<T> d1, d2;
    d1 = d;
    computed2(l, PsiLambdaCheck_vs_v[l], d, d1, d2);

    // Compute underlinePiCheck
    CoefficientsByLevel<T> dummy, PhiPiunderlineCheck_vs_v;
    trial_localtransform.reconstruct(dummy, PsiLambdaCheck_vs_v[l], l, PhiPiunderlineCheck_vs_v);
//    std::cout << "l = " << l << ", PhiPiunderlineCheck_vs_v = " << PhiPiunderlineCheck_vs_v  << std::endl;

    // Compute underline d
    CoefficientsByLevel<T> underline_d;
    trial_localtransform.reconstruct(d2, dummy, l, underline_d);
//    std::cout << "l = " << l << ", underline_d = " << underline_d << std::endl;

    //Compute <\Phi_{PiCheck1} , \Phi_{supp(underline_d)} >
    applyBilinearForm(l+1, underline_d, PhiPiunderlineCheck_vs_v);

    test_localtransform.decompose_(PhiPiunderlineCheck_vs_v, l+1,
                                   dummy, PsiLambdaCheck_vs_v[l]);

    // Compute underline d
    trial_localtransform.reconstruct(dummy, c[l], l, underline_d);

    //Recursive call of eval
    this->evalL(l+1, underline_d, c, PsiLambdaCheck_vs_v);

    return;
}

}   // namespace lawa
