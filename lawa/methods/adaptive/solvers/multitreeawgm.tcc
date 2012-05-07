namespace lawa {

template <typename Index, typename Basis, typename LocalOperator, typename RHS, typename Preconditioner>
MultiTreeAWGM<Index,Basis,LocalOperator,RHS,Preconditioner>::
MultiTreeAWGM(const Basis &_basis, LocalOperator &_A, RHS &_F, Preconditioner &_Prec)
: basis(_basis), A(_A), F(_F), Prec(_Prec), alpha(0.5), gamma(0.1), hashMapSize(SIZEHASHINDEX2D)
{

}

template <typename Index, typename Basis, typename LocalOperator, typename RHS, typename Preconditioner>
void
MultiTreeAWGM<Index,Basis,LocalOperator,RHS,Preconditioner>::setParameters(T _alpha, T _gamma)
{
    alpha = _alpha;
    gamma = _gamma;
}

template <typename Index, typename Basis, typename LocalOperator, typename RHS, typename Preconditioner>
void
MultiTreeAWGM<Index,Basis,LocalOperator,RHS,Preconditioner>::
cg_solve(Coefficients<Lexicographical,T,Index> &u, T _eps, const char *filename,
      int NumOfIterations, T EnergyNorm)
{
    Coefficients<Lexicographical,T,Index> r(hashMapSize), p(hashMapSize),
                                          Ap(hashMapSize);

    for (int iter=0; iter<=100; ++iter) {

        for (coeff_it it=r.begin(); it!=r.end(); ++it) {
            if (u.find((*it).first)==u.end()) Ap.erase((*it).first);
            else                              Ap[(*it).first] = 0.;
        }
        r.clear();
        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            if (u.find((*it).first)!=u.end()) {
                T tmp = Prec((*it).first) * F((*it).first);
                r[(*it).first] = 0.;
                p[(*it).first] = tmp;
            }
        }

        A.eval(u,r,Prec);
        r -= p;
        p = r;
        p *= (T)(-1.);

        T rNormSquare = r.norm((T)2.);
        int cg_iters=0;
        int maxIterations=100;
        T mv_time1=0., mv_time2=0.;
        for (cg_iters=0; cg_iters<maxIterations; ++cg_iters) {
            Ap.setToZero();
            A.eval(p,Ap,Prec);
            T pAp = p * Ap;
            T alpha = rNormSquare/pAp;
            p *= alpha;
            u += p;
            p *= (T)1./alpha;

            Ap *= alpha;
            r += Ap;

            T rNormSquarePrev = rNormSquare;
            rNormSquare = r.norm((T)2.);
            T beta = rNormSquare/rNormSquarePrev;
            p *= beta;
            p -= r;
        }
    }

    r.setToZero();
    extendMultiTree(basis, u, r);
    A.eval(u,r,Prec);
    for (coeff_it it=r.begin(); it!=r.end(); ++it) {
        (*it).second -= Prec((*it).first) * F((*it).first);
    }

    long double r_norm = (long double)r.norm((T)2.);
    T threshbound = std::sqrt(1-alpha*alpha) * r.norm(2.)/std::sqrt(T(r.size()));
    Coefficients<Bucket,T,Index> r_bucket;
    r_bucket.bucketsort(r, threshbound);

    long double P_Lambda_r_norm_square = 0.0L;
    if (u.size() > 0) {
        for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
            P_Lambda_r_norm_square += std::pow(r[(*it).first],(T)2.);
            r.erase((*it).first);
        }
    }


    std::cerr << "      ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_r_norm_square)
              << ", alpha*r_norm = " << alpha*r_norm << std::endl;

    for (int i=0; i<(int)r_bucket.bucket_ell2norms.size(); ++i) {
        P_Lambda_r_norm_square += std::pow(r_bucket.bucket_ell2norms[i],2.0L);
        r_bucket.addBucketToCoefficients(p,i);
        if (P_Lambda_r_norm_square >= alpha*r_norm*alpha*r_norm) break;
    }

    for (const_coeff_it it=p.begin(); it!=p.end(); ++it) {
        completeMultiTree(basis, (*it).first, u);
    }
}

}   // namespace lawa
