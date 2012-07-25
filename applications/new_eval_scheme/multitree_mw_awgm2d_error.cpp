#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef long double T;

typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

typedef Basis<T,Orthogonal,Interval,Multi>                          PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;
typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>             Basis2D;

typedef OptimizedH1Preconditioner2D<T,Basis2D>                      Preconditioner;

///  Underlying bilinear form
typedef RefinementBasis::LaplaceOperator1D                          RefinementLaplaceOp1D;
typedef AdaptiveLaplaceOperator1D<T,Orthogonal,Interval,Multi>      LaplaceOp1D;

///  Local operator in 1d
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementLaplaceOp1D,LaplaceOp1D>          LocalOp1D;

typedef UniDirectionalLocalOperator<Index2D,XOne,LocalOp1D,
                                            NotXOne,Index1D>        UniDirectionalLocalOpXOne2D;
typedef UniDirectionalLocalOperator<Index2D,XTwo,LocalOp1D,
                                            NotXTwo,Index1D>        UniDirectionalLocalOpXTwo2D;

typedef CompoundLocalOperator<Index2D, UniDirectionalLocalOpXOne2D,
                              UniDirectionalLocalOpXTwo2D>          CompoundLocalOperator2D;

//Righthandsides definitions (separable)
typedef RHSWithPeaks1D<T,PrimalBasis>                               Rhs1D;
typedef AdaptiveSeparableRhs<T,Index2D,Rhs1D,Rhs1D >                AdaptiveSeparableRhsIntegral2D;
typedef CompoundRhs<T,Index2D,AdaptiveSeparableRhsIntegral2D,
                    AdaptiveSeparableRhsIntegral2D>                 CompoundRhsIntegral2D;


typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::iterator           coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

int example = 2;

T u1(T x)   {    return 1.; }

T u2(T y)   {    return 1.; }

T du1(T x)  {    return 0.; }

T du2(T y)  {    return 0.; }

T ddu1(T x) {    return -10.;   }

T ddu2(T y) {    return -10.;   }

long double EnergyErrorSquared = 0.L;

/*
int example = 3;
T u1(T x)   {    return x*x*(1-x)*(1-x); }

T u2(T y)   {    return y*y*(1-y)*(1-y); }

T du1(T x)  {    return 2*x*(1-x)*(1-x)-2*x*x*(1-x); }

T du2(T y)  {    return 2*y*(1-y)*(1-y)-2*y*y*(1-y); }

T ddu1(T x) {    return 2*(1-x)*(1-x) - 8*x*(1-x) + 2*x*x; }

T ddu2(T y) {    return 2*(1-y)*(1-y) - 8*y*(1-y) + 2*y*y; }

long double EnergyErrorSquared = 2.*(1.L/630.L * 2.L/105.L);
*/
T f1(T x)   {   return -ddu1(x); }

T f2(T y)   {   return -ddu2(y); }

T sol(T x, T y) {   return u1(x) * u2(y); }

template <typename T>
void
setUp_f_eps(int example, PrimalBasis &basis,
            Preconditioner &Prec, Coefficients<Lexicographical,T,Index2D> &f_eps,
            Coefficients<Lexicographical,T,Index1D> &rhs_u1_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u1,
            Coefficients<Lexicographical,T,Index1D> &rhs_u2_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u2,
            Coefficients<Lexicographical,T,Index1D> &rhs_f1_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f1,
            Coefficients<Lexicographical,T,Index1D> &rhs_f2_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f2);

template <typename T>
void
test_extendMultiTree(const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D>  &v,
                     Coefficients<Lexicographical,T,Index2D>  &C_v, int coordDirection);

int main (int argc, char *argv[]) {

    cout.precision(20);
    if (argc!=4) {
        cout << "Usage: " << argv[0] << " d j0 J" << endl;
        return 0;
    }
    int d   = atoi(argv[1]);
    int j0  = atoi(argv[2]);
    int J  = atoi(argv[3]);
    T alpha = 0.7;
    T gamma = 0.1;
    const char* residualType = "standard";
    const char* treeType = "sparsetree";    // "gradedtree";
    T eps   = 1e-2;
    Timer time;

    /// Basis initialization
    PrimalBasis       basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis  &refinementbasis = basis.refinementbasis;
    Basis2D basis2d(basis,basis);

    /// Operator initialization
    LaplaceOp1D                  laplaceOp1D(basis);
    LocalOp1D                    localOp1D(basis,basis,refinementbasis.LaplaceOp1D,laplaceOp1D);
    UniDirectionalLocalOpXOne2D  uniDirectionalOpXOne2D(localOp1D);
    UniDirectionalLocalOpXTwo2D  uniDirectionalOpXTwo2D(localOp1D);
    CompoundLocalOperator2D      localOp2D(uniDirectionalOpXOne2D,uniDirectionalOpXTwo2D);

    /// Initialization of preconditioner
    Preconditioner  Prec(basis2d,1.,1.,0.);

    /// Initialization of rhs
    DenseVectorT sing_pts_x, sing_pts_y;
    DenseMatrixT no_deltas, deltas_x, deltas_y;
    int order = 20;
    if (example==2) {  int order = 4+2*d; }

    Function<T>                    fct_u1(u1,sing_pts_x), fct_f1(f1,sing_pts_x);
    Function<T>                    fct_u2(u2,sing_pts_y), fct_f2(f2,sing_pts_y);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u1(basis, fct_u1, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f1(basis, fct_f1, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u2(basis, fct_u2, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f2(basis, fct_f2, no_deltas, order);
    Coefficients<Lexicographical,T,Index1D> rhs_u1_data(SIZEHASHINDEX1D),
                                            rhs_f1_data(SIZEHASHINDEX1D),
                                            rhs_u2_data(SIZEHASHINDEX1D),
                                            rhs_f2_data(SIZEHASHINDEX1D);
    AdaptiveSeparableRhsIntegral2D rhs1(rhs_f1, rhs_f1_data, rhs_u2, rhs_u2_data);
    AdaptiveSeparableRhsIntegral2D rhs2(rhs_u1, rhs_u1_data, rhs_f2, rhs_f2_data);
    CompoundRhsIntegral2D          F(rhs1,rhs2);

    Coefficients<Lexicographical,T,Index2D> f_eps(SIZEHASHINDEX2D);
    setUp_f_eps<T>(example, basis, Prec, f_eps,
                   rhs_u1_data, rhs_u1, rhs_u2_data, rhs_u2,
                   rhs_f1_data, rhs_f1, rhs_f2_data, rhs_f2);

    //T tol = 3e-07;
    cout << "Norm of f_eps: " << f_eps.norm() << endl;
    T tol = 0.1;

    stringstream residual_error_filename;
    residual_error_filename << "error_multitree_mw_awgm_poisson2d_" << example << "_"
                            << argv[1] << "_" << argv[2] << "_" << alpha << "_" << gamma << "_"
                            << residualType << "_" << treeType<< ".dat";
    ofstream residual_error_file(residual_error_filename.str().c_str());

    Coefficients<Lexicographical,T,Index2D> u(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> r_multitree(SIZEHASHINDEX2D);

    for (int iter=1; iter<=100; ++iter) {
        Coefficients<Lexicographical,T,Index2D> r_apply(SIZEHASHINDEX2D);
        Coefficients<Lexicographical,T,Index2D> r_eps(SIZEHASHINDEX2D);
        stringstream coefffilename;
        coefffilename << "coeff2d/coeff_multitree_mw_awgm_poisson2d_" << example << "_"
                      << argv[1] << "_" << argv[2] << "_" << alpha << "_" << gamma << "_"
                      << residualType << "_" << treeType << "__" << iter << ".dat";
        readCoefficientsFromFile(u, coefffilename.str().c_str());
        if (u.size()==0) break;
        cout << "Size of u: " << u.size() << endl;
        localOp2D.apply(u,r_eps,Prec,eps);
        r_eps -= f_eps;
        cout << "Reference computation of r_eps has finished." << endl;
        T exact_residual = r_eps.norm();
        cout << "Size of r_eps = " << r_eps.size() << endl;

        time.start();
        //r = u;
        r_multitree.setToZero();
        extendMultiTree(basis2d, u, r_multitree, residualType);
        localOp2D.eval(u,r_multitree,Prec);
        for (coeff2d_it it=r_multitree.begin(); it!=r_multitree.end(); ++it) {
            (*it).second -= Prec((*it).first) * F((*it).first);
        }
        time.stop();
        cout << "Size of multitree residual = " << r_multitree.size() << endl;
        T new_residual_time = time.elapsed();
        T new_residual_norm = r_multitree.norm();
        int new_residual_length = r_multitree.size();
        //cout << "r = " << r << endl;

        r_eps -= r_multitree;
        T new_residual_diff = r_eps.norm();
        r_eps += r_multitree;
        //cout << "diff = " << r << endl;

        T apply_residual_norm = 0., apply_residual_norm1 = 0., apply_residual_norm2 = 0.;
        int apply_residual_length = 0, apply_residual_length1 = 0, apply_residual_length2 = 0;
        T apply_residual_diff = 0., apply_residual_diff1 = 0., apply_residual_diff2 = 0.;
        T apply_residual_time = 0., apply_residual_time1 = 0., apply_residual_time2 = 0.;

        while(1) {
            r_apply.setToZero();
            time.start();
            localOp2D.apply(u,r_apply,Prec,tol/2.);
            cout << "      " << r_apply.size() << endl;
            //r -= f_eps;
            r_apply -= THRESH(f_eps,tol/2.,true,true);
            time.stop();
            apply_residual_norm2 = r_apply.norm();
            apply_residual_length2 = r_apply.size();
            apply_residual_time2 = time.elapsed();
            r_eps -= r_apply;
            apply_residual_diff2 = r_eps.norm();
            r_eps += r_apply;
            cerr << "   DEBUG: tol = " << tol
                 << ", apply_residual_diff = " << apply_residual_diff2
                 << ", new_residual_diff = " << new_residual_diff
                 << ", apply_residual_time = " << time.elapsed() << endl;;
            if (apply_residual_diff2<new_residual_diff)  {
                if (apply_residual_time2 < apply_residual_time1) {
                    apply_residual_norm = apply_residual_norm2;
                    apply_residual_diff = apply_residual_diff2;
                    apply_residual_length = apply_residual_length2;
                    apply_residual_time = apply_residual_time2;
                }
                else {
                    apply_residual_norm = apply_residual_norm1;
                    apply_residual_diff = apply_residual_diff1;
                    apply_residual_length = apply_residual_length1;
                    apply_residual_time = apply_residual_time1;
                }
                break;
            }
            else {
                apply_residual_norm1 = apply_residual_norm2;
                apply_residual_diff1 = apply_residual_diff2;
                apply_residual_length1 = apply_residual_length2;
                apply_residual_time1 = apply_residual_time2;
                tol *= 0.9;
            }
        }

        cout << u.size() << " " << r_eps.norm() << ":" << endl;
        cout << "new residual:   " << new_residual_norm << " " << new_residual_diff << " " << new_residual_length << endl;
        cout << "apply residual: " << apply_residual_norm << " " << apply_residual_diff << " " << apply_residual_length << endl << endl;

        residual_error_file << u.size() << " " << exact_residual << " "
                            << new_residual_norm << " " << new_residual_diff << " "
                            << new_residual_length << " " << new_residual_time << " "
                            << apply_residual_norm << " " << apply_residual_diff << " "
                            << apply_residual_length << " " << apply_residual_time <<  endl;

        eps = min(eps, 0.1*tol);

    }
    return 0;
}

template <typename T>
void
setUp_f_eps(int example, PrimalBasis &basis,
            Preconditioner &Prec, Coefficients<Lexicographical,T,Index2D> &f_eps,
            Coefficients<Lexicographical,T,Index1D> &rhs_u1_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u1,
            Coefficients<Lexicographical,T,Index1D> &rhs_u2_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u2,
            Coefficients<Lexicographical,T,Index1D> &rhs_f1_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f1,
            Coefficients<Lexicographical,T,Index1D> &rhs_f2_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f2)
{
    int j0 = basis.j0;
    if (example==2) {
        for (int k=basis.mra.rangeIL(j0).firstIndex(); k<=basis.mra.rangeIL(j0).lastIndex(); ++k) {
            Index1D index1d(j0,k,XBSpline);
            rhs_u1_data[index1d] = rhs_u1(index1d); rhs_u2_data[index1d] = rhs_u2(index1d);
            rhs_f1_data[index1d] = rhs_f1(index1d); rhs_f2_data[index1d] = rhs_f2(index1d);
        }
        for (int k=basis.mra.rangeIR(j0).firstIndex(); k<=basis.mra.rangeIR(j0).lastIndex(); ++k) {
            Index1D index1d(j0,k,XBSpline);
            rhs_u1_data[index1d] = rhs_u1(index1d); rhs_u2_data[index1d] = rhs_u2(index1d);
            rhs_f1_data[index1d] = rhs_f1(index1d); rhs_f2_data[index1d] = rhs_f2(index1d);
        }

        for (int j=j0; j<=25; ++j) {
            for (int k=basis.rangeJL(j).firstIndex(); k<=basis.rangeJL(j).lastIndex(); ++k) {
                Index1D index1d(j,k,XWavelet);
                rhs_u1_data[index1d] = rhs_u1(index1d); rhs_u2_data[index1d] = rhs_u2(index1d);
                rhs_f1_data[index1d] = rhs_f1(index1d); rhs_f2_data[index1d] = rhs_f2(index1d);
            }
            for (int k=basis.rangeJR(j).firstIndex(); k<=basis.rangeJR(j).lastIndex(); ++k) {
                Index1D index1d(j,k,XWavelet);
                rhs_u1_data[index1d] = rhs_u1(index1d); rhs_u2_data[index1d] = rhs_u2(index1d);
                rhs_f1_data[index1d] = rhs_f1(index1d); rhs_f2_data[index1d] = rhs_f2(index1d);
            }
        }

        for (const_coeff1d_it it_x=rhs_u1_data.begin(); it_x!=rhs_u1_data.end(); ++it_x) {
            for (const_coeff1d_it it_y=rhs_u2_data.begin(); it_y!=rhs_u2_data.end(); ++it_y) {
                Index2D index((*it_x).first,(*it_y).first);
                f_eps[index] =  ( (*it_x).second * rhs_f2_data[(*it_y).first]
                               + rhs_f1_data[(*it_x).first] * (*it_y).second) * Prec(index);
            }
        }
    }
    std::cerr << "#Supp f_eps = " << f_eps.size() << std::endl;
}

/*
Coefficients<Lexicographical,T,Index2D> r1(SIZEHASHINDEX2D), r2(SIZEHASHINDEX2D);
test_extendMultiTree(basis2d, u, r1, 1);
test_extendMultiTree(basis2d, u, r2, 2);
r1 += r2;
cout << "DEBUG: size of r = " << r1.size() << endl;

template <typename T>
void
test_extendMultiTree(const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D>  &v,
                     Coefficients<Lexicographical,T,Index2D>  &C_v, int coordDirection)
{
    for (const_coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        C_v[(*it).first] = (T)0.;
        IndexSet<Index1D > C_index1, C_index2;
        if (coordDirection==1) {
            C_index1 = C((*it).first.index1, (T)1., basis.first);
            for (const_set1d_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
                Index2D newindex((*it_C_index1), (*it).first.index2);
                if (C_v.find(newindex)==C_v.end()) { completeMultiTree(basis,newindex,C_v); }
            }
        }
        else if (coordDirection==2) {
            C_index2 = C((*it).first.index2, (T)1., basis.second);
            for (const_set1d_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
                Index2D newindex((*it).first.index1, (*it_C_index2));
                if (C_v.find(newindex)==C_v.end()) { completeMultiTree(basis,newindex,C_v); }
            }
        }
        else {
            cerr << "Non-admissible coordinate direction." << endl;
            exit(1);
        }
    }
}
*/

/*
        Index2D maxIndex;
        Index2D maxWaveletIndex;
        int *jmax = new int[2];
        int arrayLength = 2;
        getLevelInfo(u, maxIndex, maxWaveletIndex, jmax, arrayLength);
        int J = max(jmax[0],jmax[1]);
        extendMultiTreeAtBoundary(basis2d, u, r, J+3);
        */
