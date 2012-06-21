namespace lawa {

template <typename Index, typename Basis, typename LocalOperator, typename RHS, typename Preconditioner>
MultiTreeAWGM<Index,Basis,LocalOperator,RHS,Preconditioner>::
MultiTreeAWGM(const Basis &_basis, LocalOperator &_A, RHS &_F, Preconditioner &_Prec,
              Coefficients<Lexicographical,T,Index> &_f_eps)
: basis(_basis), A(_A), F(_F), Prec(_Prec), f_eps(_f_eps),
  alpha(0.5), gamma(0.1), residualType("standard"), hashMapSize(SIZEHASHINDEX2D),
  compute_f_minus_Au_error(false), write_coefficients_to_file(false)
{

}

template <typename Index, typename Basis, typename LocalOperator, typename RHS, typename Preconditioner>
void
MultiTreeAWGM<Index,Basis,LocalOperator,RHS,Preconditioner>::setParameters(T _alpha, T _gamma,
                                                                           const char* _residualType,
                                                                           bool _compute_f_minus_Au_error,
                                                                           bool _write_coefficients_to_file)
{
    alpha = _alpha;
    gamma = _gamma;
    residualType = _residualType;
    compute_f_minus_Au_error = _compute_f_minus_Au_error;
    write_coefficients_to_file = _write_coefficients_to_file;
    std::cerr << "Residual Type: " << residualType << std::endl;
}

template <typename Index, typename Basis, typename LocalOperator, typename RHS, typename Preconditioner>
void
MultiTreeAWGM<Index,Basis,LocalOperator,RHS,Preconditioner>::
cg_solve(Coefficients<Lexicographical,T,Index> &u, T _eps, int NumOfIterations,
         T EnergyNormSquared, const char *filename, const char *coefffilename)
{
    Coefficients<Lexicographical,T,Index> r(hashMapSize),       // residual vector for cg
                                          p(hashMapSize),       // auxiliary vector for cg
                                          Ap(hashMapSize);      // auxiliary vector for cg
    Coefficients<Lexicographical,T,Index> res(hashMapSize);     // approximate residual for f-Au
    Coefficients<Lexicographical,T,Index> u_leafs(hashMapSize); // "leafs" of u

    Index maxIndex;
    Index maxWaveletIndex;
    int *jmax = new int[1];
    int arrayLength = 1;
    long double Residual = 1.L;

    Timer time;
    Timer iteration_time;

    for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
        u_leafs[(*it).first] = 0.;
        r[(*it).first] = 0.;
        p[(*it).first] = 0.;
        Ap[(*it).first] = 0.;
        res[(*it).first] = 0.;
    }

    std::ofstream file(filename);
    file.precision(16);

    T time_total_comp = 0.;
    for (int iter=1; iter<=NumOfIterations; ++iter) {
        T time_mv_linsys = 0.;
        T time_mv_residual = 0.;
        T time_multitree_residual = 0.;
        iteration_time.start();

        int N = u.size();
        int N_residual = 0;


        std::cerr << std::endl << "*****  Iteration " << iter << " *****" << std::endl << std::endl;
        std::cerr << "   Current number of dof = " << u.size() << std::endl;

        /* ******************* Resetting of vectors *********************** */

        std::cerr << "   Resetting vectors..." << std::endl;

        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            T tmp = Prec((*it).first) * F((*it).first);
            r[(*it).first] = 0.;
            p[(*it).first] = tmp;
        }
        std::cerr.precision(16);
        std::cerr << "   || f ||_{ell_2} = " << p.norm((T)2.) << std::endl;
        std::cerr.precision(6);
        std::cerr << "   ...finished." << std::endl;


        /* ******************* CG method for Galerkin system *********************** */

        std::cerr << "   CG method started." << std::endl;
        A.eval(u,r,Prec,"galerkin");
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
            time.start();
            A.eval(p,Ap,Prec,"galerkin");
            time.stop();
            std::cerr << "      DEBUG: " << time.elapsed() << std::endl;
            time_mv_linsys += time.elapsed();
            T pAp = p * Ap;
            T alpha = cg_rNormSquare/pAp;
            p *= alpha;
            u += p;
            p *= (T)1./alpha;
            Ap *= alpha;
            r += Ap;

            T cg_rNormSquarePrev = cg_rNormSquare;
            cg_rNormSquare = r*r;
            std::cerr << "      Current error in cg: " << std::sqrt(cg_rNormSquare) << std::endl;
            T beta = cg_rNormSquare/cg_rNormSquarePrev;
            p *= beta;
            p -= r;
        }
        std::cerr << "   CG method finished after " << cg_iter << " iterations." << std::endl;
        time_mv_linsys *= 1./cg_iter;
        if (write_coefficients_to_file) {
            writeCoefficientsToFile(u, iter, coefffilename);
        }

        /* ************************************************************************* */

        /* ***************** Errors and some info on the solution ******************** */

        std::cerr << "   Computing energy error..." << std::endl;
        Ap.setToZero();
        A.eval(u,Ap,Prec,"galerkin");
        T uAu = Ap*u;
        T fu  = 0.;
        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            fu += (*it).second*Prec((*it).first) * F((*it).first);
        }
        T EnergyError =  sqrt(fabs(EnergyNormSquared - 2*fu + uAu));
        std::cerr << "   ... finished with energy error: " << EnergyError << std::endl;

        T f_minus_Au_error = 0.;
        if (compute_f_minus_Au_error) {
            compute_f_minus_Au(u, 1e-9, f_minus_Au_error);
        }


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


        /* ************************************************************************* */

        /* ******************* Computing approximate residual ********************** */

        std::cerr << "   Computing residual..." << std::endl;
        res.setToZero();
        std::cerr << "     Computing multi-tree for residual..." << std::endl;
        std::cerr << "        #supp u = " << u.size() << ", #supp r = " << res.size() << std::endl;
        time.start();
        extendMultiTree(basis, u_leafs, res, residualType);
        //extendMultiTreeAtBoundary(basis, u, res, J+3);
        time.stop();
        time_multitree_residual = time.elapsed();
        N_residual = res.size();
        std::cerr << "     ... finished after " << time.elapsed() << std::endl;
        std::cerr << "   #supp u = " << u.size() << ", #supp r = " << res.size() << std::endl;
        std::cerr << "     Computing matrix vector product..." << std::endl;
        time.start();
        A.eval(u,res,Prec);
        std::cerr << "     ... finished." << std::endl;
        std::cerr << "     Substracting right-hand side..." << std::endl;
        for (coeff_it it=res.begin(); it!=res.end(); ++it) {
            (*it).second -= Prec((*it).first) * F((*it).first);
        }
        time.stop();
        time_mv_residual = time.elapsed();
        time_multitree_residual += time_mv_residual;
        Residual = res.norm(2.);
        std::cerr << "   ... finished with Residual: " << Residual << std::endl;

        /* ************************************************************************* */

        /* ********************** Computing next index set  ************************ */

        long double P_Lambda_Residual_square = 0.0L;
        if (u.size() > 0) {
            for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
                P_Lambda_Residual_square += std::pow(r[(*it).first],(T)2.);
                res.erase((*it).first);
            }
        }
        if (res.size()!=0) {
            T threshbound = std::sqrt(1-alpha*alpha) * res.norm((T)2.)/std::sqrt(T(res.size()));
            Coefficients<Bucket,T,Index> r_bucket;
            r_bucket.bucketsort(res, threshbound);
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

        // Above we set res = res-res|_{supp u}. Now we change this back.
        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            res[(*it).first] = 0.;
        }

        // The vector p satisfies the bulk criterion. But it is not yet a multi-tree...
        for (const_coeff_it it=p.begin(); it!=p.end(); ++it) {
            completeMultiTree(basis, (*it).first, u);
        }

        // Note that r has still supported on the previous Galerkin index set!
        u_leafs.clear();
        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            if (r.find((*it).first)==r.end()) u_leafs[(*it).first] = 0.;
        }

        /* ************************************************************************* */

        iteration_time.stop();
        time_total_comp += iteration_time.elapsed();

        file << N << " " << N_residual << " " << time_total_comp << " " << EnergyError
                  << " " << Residual << " " << f_minus_Au_error << " "
                  << time_mv_linsys << " " << time_mv_residual << " " << time_multitree_residual << std::endl;

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
    filenamestr << filename << "__" << i << ".dat";
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
