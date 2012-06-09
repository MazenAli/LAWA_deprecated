namespace lawa {

template <typename Index, typename Basis, typename LocalOperator, typename RHS, typename Preconditioner>
MultiTreeAWGM<Index,Basis,LocalOperator,RHS,Preconditioner>::
MultiTreeAWGM(const Basis &_basis, LocalOperator &_A, RHS &_F, Preconditioner &_Prec,
              Coefficients<Lexicographical,T,Index> &_f_eps)
: basis(_basis), A(_A), F(_F), Prec(_Prec), f_eps(_f_eps),
  alpha(0.5), gamma(0.1), residualType("standard"), hashMapSize(SIZEHASHINDEX2D),
  compute_f_minus_Au_error(false)
{

}

template <typename Index, typename Basis, typename LocalOperator, typename RHS, typename Preconditioner>
void
MultiTreeAWGM<Index,Basis,LocalOperator,RHS,Preconditioner>::setParameters(T _alpha, T _gamma,
                                                                           const char* _residualType,
                                                                           bool _compute_f_minus_Au_error)
{
    alpha = _alpha;
    gamma = _gamma;
    residualType = _residualType;
    compute_f_minus_Au_error = _compute_f_minus_Au_error;
    std::cerr << "Residual Type: " << residualType << std::endl;
}

template <typename Index, typename Basis, typename LocalOperator, typename RHS, typename Preconditioner>
void
MultiTreeAWGM<Index,Basis,LocalOperator,RHS,Preconditioner>::
cg_solve(Coefficients<Lexicographical,T,Index> &u, T _eps, const char *filename,
      int NumOfIterations, T EnergyNormSquared)
{
    Coefficients<Lexicographical,T,Index> r(hashMapSize), p(hashMapSize),
                                          Ap(hashMapSize);

    Index maxIndex;
    Index maxWaveletIndex;
    int *jmax = new int[1];
    int arrayLength = 1;
    long double Residual = 1.L;

    Timer time;
    Timer iteration_time;

    for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
        r[(*it).first] = 0.;
        p[(*it).first] = 0.;
        Ap[(*it).first] = 0.;
    }

    std::ofstream file(filename);
    file.precision(16);

    T time_total_comp = 0.;
    for (int iter=1; iter<=NumOfIterations; ++iter) {
        T time_mv_linsys = 0.;
        T time_mv_residual = 0.;
        iteration_time.start();

        int N = u.size();
        int N_residual = 0;


        std::cerr << std::endl << "*****  Iteration " << iter << " *****" << std::endl << std::endl;
        std::cerr << "   Current number of dof = " << u.size() << std::endl;
        std::cerr << "   Resetting vectors..." << std::endl;
        for (coeff_it it=r.begin(); it!=r.end(); ++it) {
            if (u.find((*it).first)==u.end()) Ap.erase((*it).first);
            else                              Ap[(*it).first] = 0.;
        }
        r.clear();
        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            T tmp = Prec((*it).first) * F((*it).first);
            r[(*it).first] = 0.;
            p[(*it).first] = tmp;
        }
        std::cerr.precision(16);
        std::cerr << "   || f ||_{ell_2} = " << p.norm((T)2.) << std::endl;
        std::cerr.precision(6);
        std::cerr << "   ...finished." << std::endl;

        std::cerr << "   CG method started." << std::endl;
        A.eval(u,r,Prec);
        r -= p;
        p = r;
        p *= (T)(-1.);
        T cg_rNormSquare = r*r;
        T tol = std::min((T)1e-2,gamma*(T)Residual);
        int maxIterations=100;
        int cg_iter=0;
        for (cg_iter=0; cg_iter<maxIterations; ++cg_iter) {
            if (std::sqrt(cg_rNormSquare)<=tol) {
                std::cerr << "      CG stopped with error " << sqrt(cg_rNormSquare) << std::endl;
                break;
            }
            std::cerr << "    Iteration " << cg_iter+1 << std::endl;
            Ap.setToZero();
            //std::cerr << "p = " << p << std::endl;
            time.start();
            A.eval(p,Ap,Prec);
            time.stop();
            time_mv_linsys += time.elapsed();
            //std::cerr << "Ap = " << Ap << std::endl;
            T pAp = p * Ap;
            T alpha = cg_rNormSquare/pAp;
            //std::cerr << "alpha = " << alpha << std::endl;
            p *= alpha;
            //std::cerr << "alpha*p = " << p << std::endl;
            u += p;
            //std::cerr << "u = " << u << std::endl;
            p *= (T)1./alpha;

            Ap *= alpha;
            r += Ap;
            //std::cerr << "r = " << r << std::endl;

            T cg_rNormSquarePrev = cg_rNormSquare;
            cg_rNormSquare = r*r;
            std::cerr << "      Current error in cg: " << std::sqrt(cg_rNormSquare) << std::endl;
            T beta = cg_rNormSquare/cg_rNormSquarePrev;
            p *= beta;
            p -= r;
        }

        std::cerr << "   CG method finished after " << cg_iter << " iterations." << std::endl;

        time_mv_linsys *= 1./cg_iter;

        writeCoefficientsToFile(u, iter, "u");

        std::cerr << "   Computing enery error..." << std::endl;
        Ap.setToZero();
        A.eval(u,Ap,Prec);
        T uAu = Ap*u;
        T fu  = 0.;
        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            fu += (*it).second*Prec((*it).first) * F((*it).first);
        }
        T EnergyError =  sqrt(fabs(EnergyNormSquared - 2*fu + uAu));
        std::cerr << "   ... finished with EnergyError: " << EnergyError << std::endl;

        getLevelInfo(u, maxIndex, maxWaveletIndex, jmax, arrayLength);
        int J = -100;
        for (int i=0; i<arrayLength; ++i) {
            if (J<jmax[i]) J=jmax[i];
        }

        std::cerr << "   Level information: " << std::endl;
        std::cerr << "      maxIndex = " << maxIndex << ", maxWaveletIndex = "
                 << maxWaveletIndex << std::endl;
        std::cerr << "      arrayLength = " << arrayLength << std::endl;
        std::cerr << "      Highest level per coordinate direction:";
        for (int i=0; i<arrayLength; ++i) {
            std::cerr << " " << jmax[i];
        }
        std::cerr << std::endl;

        //bool useSupportCenter=true;
        //plotScatterCoeff(u, basis, "u_coeff", useSupportCenter);
        //std::cerr << "Please hit enter." << std::endl;
        //getchar();

        std::cerr << "   Computing residual..." << std::endl;
        r.setToZero();
        std::cerr << "     Computing multi-tree for residual..." << std::endl;
        //Coefficients<Lexicographical,T,Index> tmp(hashMapSize);
        extendMultiTree(basis, u, r, residualType);
        //extendMultiTree(basis, tmp, r);
        //tmp.clear();
        //extendMultiTreeAtBoundary(basis, u, r, J+1);
        N_residual = r.size();
        std::cerr << "     ... finished." << std::endl;
        std::cerr << "   #supp u = " << u.size() << ", #supp r = " << r.size() << std::endl;
        std::cerr << "     Computing matrix vector product..." << std::endl;
        time.start();
        A.eval(u,r,Prec);
        time.stop();
        time_mv_residual = time.elapsed();
        std::cerr << "     ... finished." << std::endl;
        std::cerr << "     Substracting right-hand side..." << std::endl;
        time.start();
        for (coeff_it it=r.begin(); it!=r.end(); ++it) {
            (*it).second -= Prec((*it).first) * F((*it).first);
        }
        Residual = r.norm(2.);
        time.stop();
        std::cerr << "   ... finished with Residual: " << Residual << std::endl;


        T f_minus_Au_error = 0.;
        std::cerr << "#Supp u = " << u.size() << std::endl;
        if (compute_f_minus_Au_error) { compute_f_minus_Au(u, 1e-9, f_minus_Au_error); }
        std::cerr << "#Supp u = " << u.size() << std::endl;

        long double P_Lambda_Residual_square = 0.0L;
        if (u.size() > 0) {
            for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
                P_Lambda_Residual_square += std::pow(r[(*it).first],(T)2.);
                r.erase((*it).first);
            }
        }
        if (r.size()!=0) {
            T threshbound = std::sqrt(1-alpha*alpha) * r.norm((T)2.)/std::sqrt(T(r.size()));
            Coefficients<Bucket,T,Index> r_bucket;
            r_bucket.bucketsort(r, threshbound);
            std::cerr << "      ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_Residual_square)
                      << ", alpha*Residual = " << alpha*Residual << std::endl;

            for (int i=0; i<(int)r_bucket.bucket_ell2norms.size(); ++i) {
                P_Lambda_Residual_square += std::pow(r_bucket.bucket_ell2norms[i],2.0L);
                r_bucket.addBucketToCoefficients(p,i);
                if (P_Lambda_Residual_square >= alpha*Residual*alpha*Residual) {
                    //r_bucket.addBucketToCoefficients(p,i+1);
                    break;
                }
            }
        }

        for (const_coeff_it it=p.begin(); it!=p.end(); ++it) {
            completeMultiTree(basis, (*it).first, u);
        }
        iteration_time.stop();
        time_total_comp += iteration_time.elapsed();
        file << N << " " << N_residual << " " << time_total_comp << " " << EnergyError
                  << " " << Residual << " " << f_minus_Au_error << " "
                  << time_mv_linsys << " " << time_mv_residual << std::endl;

    }
    file.close();
}

template <typename Index, typename Basis, typename LocalOperator, typename RHS, typename Preconditioner>
void
MultiTreeAWGM<Index,Basis,LocalOperator,RHS,Preconditioner>::
compute_f_minus_Au(Coefficients<Lexicographical,T,Index> &u,
                   /*Coefficients<Lexicographical,T,Index> &Au,*/ T eps, T &f_minus_Au_error)
{
    Coefficients<Lexicographical,T,Index> Au(hashMapSize);
    std::cerr << "      APPLY started..." << std::endl;
    A.apply(u,Au,Prec,eps);
    std::cerr << "      ... finished." << std::endl;
    Au -= f_eps;
    f_minus_Au_error = Au.norm(2.);
    //Au.clear();
    std::cerr << "DEBUG: #buckets in Ap = " << Au.bucket_count() << std::endl;
    return;// Au.norm(2.);
}

template <typename Index, typename Basis, typename LocalOperator, typename RHS, typename Preconditioner>
void
MultiTreeAWGM<Index,Basis,LocalOperator,RHS,Preconditioner>::
writeCoefficientsToFile(Coefficients<Lexicographical,T,Index> &u, int i, const char* filename)
{
    std::stringstream filenamestr;
    filenamestr << filename << "_" << i << ".dat";
    std::ofstream file(filenamestr.str().c_str());
    file.precision(20);
    std::cerr << "Started writing into file..." << std::endl;
    for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
        file << (*it).first << " " << (*it).second << std::endl;
    }
    file.close();
    std::cerr << "... finished." << std::endl;
}

}   // namespace lawa
