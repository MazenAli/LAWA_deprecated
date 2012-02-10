#include <iostream>
#include <lawa/lawa.h>
#include <applications/new_eval_scheme/source/loc_single_scale_transforms.h>
#include <applications/new_eval_scheme/source/new_eval.h>
#include <applications/new_eval_scheme/source/localoperator.h>
#include <applications/new_eval_scheme/source/localoperator2d.h>
#include <applications/new_eval_scheme/source/multitreeoperations.h>
#include <lawa/methods/adaptive/datastructures/alignedindexset.h>
#include <lawa/methods/adaptive/datastructures/alignedcoefficients.h>
#include <lawa/methods/adaptive/datastructures/treecoefficients1d.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >   DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                         DenseVectorT;

typedef Basis<T,Primal,Interval,Dijkema>                            PrimalBasis;
typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>             Basis2D;

typedef AdaptiveLaplaceOperator1D<T,Primal,Interval,Dijkema>        LaplaceBilinearForm;
typedef AdaptiveIdentityOperator1D<T,Primal,Interval,Dijkema>       IdentityBilinearForm;
typedef HelmholtzOperator2D<T,Basis2D>                              HelmholtzBilinearForm2D;
typedef DiagonalMatrixPreconditioner2D<T,Basis2D,
                                       HelmholtzBilinearForm2D>     Preconditioner;

typedef LocalOperator<PrimalBasis,PrimalBasis,LaplaceBilinearForm,
                      Preconditioner>                               LocalLaplaceOp1D;
typedef LocalOperator<PrimalBasis,PrimalBasis,IdentityBilinearForm,
                      Preconditioner>                               LocalIdentityOp1D;
typedef LocalOperator2D<PrimalBasis,LocalLaplaceOp1D,
                        LocalIdentityOp1D>                          LocalLaplaceIdentityOp2D;
typedef LocalOperator2D<PrimalBasis,LocalIdentityOp1D,
                        LocalLaplaceOp1D>                           LocalIdentityLaplaceOp2D;

//Righthandsides definitions (separable)
typedef SeparableRHS2D<T,Basis2D >                                  SeparableRhsIntegral2D;

typedef SumOfTwoRHSIntegrals<T,Index2D,SeparableRhsIntegral2D,
                             SeparableRhsIntegral2D>                SumOfSeparableRhsIntegral2D;

typedef RHS<T,Index2D,SumOfSeparableRhsIntegral2D,
            Preconditioner>                                         SumOfSeparableRhs;

typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

void
writeIndexSetToFile(const IndexSet<Index2D> &Lambda, int example, int d, T threshTol, int ell, int nr);

void
readIndexSetFromFile(IndexSet<Index2D> &Lambda,  int example, int d, T threshTol, int ell, int nr);

void
getSparseGridIndexSet(const PrimalBasis &basis, IndexSet<Index2D> &Lambda, int j, T gamma=0.);

void
extendRHSIndexSet(const PrimalBasis &basis, IndexSet<Index2D> &Lambda, int j);

void
mv(LocalLaplaceIdentityOp2D &localLaplaceIdentityOp2D,
   LocalIdentityLaplaceOp2D &localIdentityLaplaceOp2D,
   Coefficients<Lexicographical,T,Index2D> &P,
   const IndexSet<Index2D> &Lambda, Coefficients<Lexicographical,T,Index2D> &v,
   Coefficients<Lexicographical,T,Index2D> &intermediate,
   Coefficients<Lexicographical,T,Index2D> &LIIAv, Coefficients<Lexicographical,T,Index2D> &IAUIv,
   Coefficients<Lexicographical,T,Index2D> &Av, T &time1, T &time2);


const T c = 0.;

/*
int example = 1;
T u1(T x)
{
    if      (0<=x && x<1./3.)    return 10*exp(x)-10.;
    else if (1./3<x && x<=2./3.) return -10.+10.*exp(1./3.)+(-10.+(30.+90.*(-(2./3.)+x))*(-(1./3.)+x))*(-(1./3.)+x);
    else                         return 10*exp(-(x-1.))-10.;
}

T u2(T y)
{
    if      (0<=y && y<1./3.)    return 10*exp(y)-10.;
    else if (1./3<y && y<=2./3.) return -10.+10.*exp(1./3.)+(-10.+(30.+90.*(-(2./3.)+y))*(-(1./3.)+y))*(-(1./3.)+y);
    else                         return 10*exp(-(y-1.))-10.;
}

T du1(T x) {
    if      (0<=x && x<1./3.)   return 10*exp(x);
    else if (1./3<x && x<2./3.) return 20.-180.*x + 270*x*x;
    else                        return -10*exp(-(x-1.));
}

T du2(T y) {
    if      (0<=y && y<1./3.)   return 10*exp(y);
    else if (1./3<y && y<2./3.) return 20.-180.*y + 270*y*y;
    else                        return -10*exp(-(y-1.));
}

T ddu1(T x)
{
    if      (0<=x && x<1./3.)   return 10*exp(x);
    else if (1./3<x && x<2./3.) return -180. + 540*x;
    else                        return 10*exp(-(x-1.));
}

T ddu2(T y)
{
    if      (0<=y && y<1./3.)   return 10*exp(y);
    else if (1./3<y && y<2./3.) return -180. + 540*y;
    else                        return 10*exp(-(y-1.));
}
T L2norm_x_sq = 10./567.*(23011.-26775.*exp(1./3.)+7560.*exp(2./3.));
T L2norm_y_sq = L2norm_x_sq;
T H1semi_x_sq = 20./3.*(-11.+15.*exp(2./3.));
T H1semi_y_sq = H1semi_x_sq;
T H1seminorm_squared = L2norm_x_sq*H1semi_y_sq + H1semi_x_sq*L2norm_y_sq;
*/

int example = 2;
T u1(T x)
{
    return 1.;
}

T u2(T y)
{
    return 1.;
}

T ddu1(T x)
{
    return -10.;
}

T ddu2(T y)
{
    return -10.;
}

T H1seminorm_squared = 14.05770149526849;


T f1(T x) {
    return -ddu1(x) + 0.5*c* u1(x);
}

T f2(T y) {
    return -ddu2(y) + 0.5*c* u2(y);
}

T sol(T x, T y)
{
    return u1(x) * u2(y);
}

int main (int argc, char *argv[]) {


    cout.precision(20);
    if (argc!=5) {
        cout << "Usage: " << argv[0] << " d d_ j0 J" << endl;
        return 0;
    }
    int d   = atoi(argv[1]);
    int d_  = atoi(argv[2]);
    int j0  = atoi(argv[3]);
    int J  = atoi(argv[4]);
    bool withDirichletBC=true;
    bool adaptive=true;
    T threshTol = 0.4;
    T r_norm = 0.1;
    T gamma = 0.2;
    int ell=1;
    Timer time;

    PrimalBasis       basis(d,d_,j0);
    if (withDirichletBC)    basis.enforceBoundaryCondition<DirichletBC>();
    Basis2D basis2d(basis,basis);

    LaplaceBilinearForm      LaplaceBil(basis);
    IdentityBilinearForm     IdentityBil(basis);
    HelmholtzBilinearForm2D  HelmholtzBil2D(basis2d,0.);
    Preconditioner           Prec(HelmholtzBil2D);
    int offset=5;
    if (d==2 && d_==2) {
        offset=2;
    }

    DenseVectorT sing_pts_x, sing_pts_y;
    DenseMatrixT no_deltas, deltas_x, deltas_y;
    if (example==1) {
        sing_pts_x.engine().resize(2); sing_pts_x(1) = 1./3.; sing_pts_x(2) = 2./3.;
        sing_pts_y.engine().resize(2); sing_pts_y(1) = 1./3.; sing_pts_y(2) = 2./3.;
        deltas_x.engine().resize(2,2); deltas_x(1,1) = 1./3.; deltas_x(1,2) = 10.+10.*exp(1./3.);
                                       deltas_x(2,1) = 2./3.; deltas_x(2,2) = 20.+10.*exp(1./3.);
        deltas_y.engine().resize(2,2); deltas_y(1,1) = 1./3.; deltas_y(1,2) = 10.+10.*exp(1./3.);
                                       deltas_y(2,1) = 2./3.; deltas_y(2,2) = 20.+10.*exp(1./3.);
    }

    SeparableFunction2D<T> SepFunc1(f1, sing_pts_x, u2, sing_pts_y);
    SeparableFunction2D<T> SepFunc2(u1, sing_pts_x, f2, sing_pts_y);
    int order = 40;
    if (example==2) order=3*d;

    SeparableRhsIntegral2D      rhsintegral_x(basis2d, SepFunc1, deltas_x, no_deltas, order);
    SeparableRhsIntegral2D      rhsintegral_y(basis2d, SepFunc2, no_deltas, deltas_y, order);
    SumOfSeparableRhsIntegral2D rhsintegral2d(rhsintegral_x,rhsintegral_y);
    SumOfSeparableRhs           F(rhsintegral2d,Prec);

    ofstream file("sparsegridconv.txt");
    file.precision(16);

    Coefficients<Lexicographical,T,Index2D> u(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> f(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> r(SIZEHASHINDEX2D),
                                            p(SIZEHASHINDEX2D),
                                            Ap(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> P(SIZEHASHINDEX2D);

    Coefficients<Lexicographical,T,Index2D> intermediate(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> LIIAv(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> IAUIv(SIZEHASHINDEX2D);

    IndexSet<Index2D> Lambda;
    getSparseGridIndexSet(basis,Lambda,0,gamma);
    //readIndexSetFromFile(Lambda,example,d,threshTol,1,13);

    LocalLaplaceOp1D            localLaplaceOp1D(basis, withDirichletBC, basis, withDirichletBC, offset, LaplaceBil, Prec, 1);
    LocalIdentityOp1D           localIdentityOp1D(basis,  withDirichletBC, basis, withDirichletBC, offset, IdentityBil,  Prec, 0);
    LocalLaplaceIdentityOp2D    localLaplaceIdentityOp2D(basis,localLaplaceOp1D,localIdentityOp1D);
    LocalIdentityLaplaceOp2D    localIdentityLaplaceOp2D(basis,localIdentityOp1D,localLaplaceOp1D);
    localLaplaceIdentityOp2D.setJ(9);
    localIdentityLaplaceOp2D.setJ(9);


    for (int iter=0; iter<=25; ++iter) {

        //readIndexSetFromFile(Lambda,example,d,threshTol,1,iter);
        //writeIndexSetToFile(Lambda,example,d,threshTol,ell,iter);
        //Lambda.clear();
        //getSparseGridIndexSet(basis,Lambda,iter,gamma);

        if (Lambda.size()>400000) break;
        cout << endl;

        cout << "******** Iteration " << iter << endl;

        cout << "   Current size of Lambda " << Lambda.size() << endl;
        f = F(Lambda);

        IndexSet<Index1D> Lambda_x, Lambda_y;
        split(Lambda,Lambda_x,Lambda_y);
        int jmin_x=0, jmax_x=0, jmin_y=0, jmax_y=0;
        getMinAndMaxLevel(Lambda_x,jmin_x,jmax_x);
        getMinAndMaxLevel(Lambda_y,jmin_y,jmax_y);
        cout << "   Max level = (" << jmax_x << ", " << jmax_y << ")" << endl;

        r.clear();
        p.clear();
        Ap.clear();
        FillWithZeros(Lambda,r);
        FillWithZeros(Lambda,p);
        FillWithZeros(Lambda,Ap);

        int maxIterations = 200;
        T tol;
        if (adaptive) tol=std::min(1e-2,1e-2*r_norm);
        else          tol=1e-8;

        T alpha, beta, rNormSquare, rNormSquarePrev;
        T dummy1=0., dummy2=0.;
        cerr << "   Computing preconditioner." << endl;
        for (const_set2d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
            if (P.find((*it))==P.end()) P[(*it)] = Prec(*it);
        }
        cerr << "   Computing matrix vector product for initial residual." << endl;

        mv(localLaplaceIdentityOp2D, localIdentityLaplaceOp2D, P,
           Lambda, u, intermediate, LIIAv, IAUIv, r, dummy1, dummy2);

        r -= f;
        p -= r;
        rNormSquare = r*r;

        int cg_iters=0;
        T mv_time1=0., mv_time2=0.;
        for (cg_iters=0; cg_iters<maxIterations; ++cg_iters) {
            if (sqrt(rNormSquare)<=tol) {
                cerr << "      CG stopped with error " << sqrt(rNormSquare) << endl;
                break;
            }
            T time1=0., time2=0.;
            time.start();
            cerr << "   Computing matrix vector product for initial residual." << endl;
            mv(localLaplaceIdentityOp2D, localIdentityLaplaceOp2D, P,
                       Lambda, p, intermediate, LIIAv, IAUIv, Ap, time1, time2);
            time.stop();
            mv_time1 += time1;
            mv_time2 += time2;
            T time_mv = time.elapsed();

            T pAp = p * Ap;
            alpha = rNormSquare/pAp;
            u += alpha*p;
            r += alpha*Ap;

            rNormSquarePrev = rNormSquare;
            rNormSquare = r*r;
            time.stop();
            cerr << "      Current error in cg: " << std::sqrt(rNormSquare) << endl;
            beta = rNormSquare/rNormSquarePrev;
            p *= beta;
            p -= r;
        }

        mv(localLaplaceIdentityOp2D, localIdentityLaplaceOp2D, P,
           Lambda, u, intermediate, LIIAv, IAUIv, Ap, dummy1, dummy2);
        T uAu = Ap*u;
        T fu  = f*u;
        T H1error =  sqrt(fabs(H1seminorm_squared - 2*fu + uAu));
        //T H1error =  sqrt(fabs(H1seminorm_squared - fu));
        std::cerr << "---> H1-error: " << H1error << " " << fu << endl;

        /*
        stringstream coeff_filename;
        if (adaptive) {
            coeff_filename << "u_adap_" << example << "_" << d << "_" << threshTol << "_" << ell << "_" << u.size();
        }
        else {
            coeff_filename << "u_sg_" << example << "_" << d << "_" << gamma << "_" << ell << "_" << u.size();
        }
        plotScatterCoeff2D<T,Index2D,PrimalBasis,PrimalBasis>(u, basis, basis, coeff_filename.str().c_str());
        */

        IndexSet<Index2D> checkLambda = Lambda;
        IndexSet<Index2D> C_Lambda = Lambda;
        for (int i=1; i<=ell; ++i) {
            C_Lambda = C(C_Lambda, 1., basis2d, true);
            for (const_set2d_it it=C_Lambda.begin(); it!=C_Lambda.end(); ++it) {
                extendMultiTree2(basis2d,(*it),offset,checkLambda);
            }
        }
        cerr << "   Residual level l = " << ell
             << ", size of of checkLambda = " << checkLambda.size() << endl;

        for (const_set2d_it it=checkLambda.begin(); it!=checkLambda.end(); ++it) {
            if (P.find((*it))==P.end()) P[(*it)] = Prec(*it);
        }

        if (example==1) {
            f = F(checkLambda);
        }
        else {
            IndexSet<Index2D> LambdaBoundary;
            extendRHSIndexSet(basis, LambdaBoundary, std::max(jmax_x,jmax_y)+2);
            f = F(LambdaBoundary);
        }
        cerr << "   Extension of rhs finished." << endl;
        T time_res1=0., time_res2=0.;
        mv(localLaplaceIdentityOp2D, localIdentityLaplaceOp2D, P,
           checkLambda, u, intermediate, LIIAv, IAUIv, r, time_res1, time_res2);
        cerr << "   Matrix vector for residual computation finished." << endl;
        r-=f;
        r_norm = r.norm(2.);
        file << Lambda.size() << " " << checkLambda.size()
             << " " << fu << " " << r_norm
             << " " << mv_time1/cg_iters << " " << mv_time2/cg_iters << " " << time_res1 << " " << time_res2 << endl;
        cerr << "---> H1-error = " << H1error << ", res = " << r_norm << endl;

        if (adaptive) {
            T P_Lambda_r_norm_square = 0.;
            for (const_set2d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
                P_Lambda_r_norm_square += std::pow(r[(*it)],2.);
                r.erase((*it));
            }
            r = THRESH(r,threshTol*r.norm(2.));
            for (const_coeff2d_it it=r.begin(); it!=r.end(); ++it) {
                extendMultiTree2(basis2d,(*it).first,offset,Lambda);
            }
        }
        else {
            getSparseGridIndexSet(basis,Lambda,iter+1,gamma);
        }
        cerr << "   Computation of new Lambda finished." << endl;
        r.clear();
    }
    std::cerr << "Start plot with u.size() == " << u.size() << std::endl;
    plot2D<T,Basis2D,Preconditioner>(basis2d, u, Prec, sol, 0., 1., 0., 1., 0.01, "sol");
    std::cerr << "Finished plot with u.size() == " << u.size() << std::endl;
    return 0;
}

void
writeIndexSetToFile(const IndexSet<Index2D> &Lambda, int example, int d, T threshTol, int ell, int nr)
{
    stringstream filename;
    filename << "_Lambda2d_" << example << "_" << d << "_" << threshTol << "_" << ell << "_" << nr << ".dat";
    ofstream file(filename.str().c_str());
    for (const_set2d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
        file << *it << endl;
    }
    file.close();
}

void
readIndexSetFromFile(IndexSet<Index2D> &Lambda, int example, int d, T threshTol, int ell, int nr)
{
    stringstream filename;
    filename << "indexsets/Lambda2d_" << example << "_" << d << "_"
             << threshTol << "_" << ell << "_" << nr << ".dat";
    std::ifstream infile (filename.str().c_str());
    if (infile.is_open()) {
        cerr << "   Indexset file is open." << endl;
    }
    else {
        cerr << "   Indexset file " << filename.str().c_str()  << " is not open." << endl;
    }

    std::string line;
    std::string field1, field2, field3, field4, field5, field6;
    while(std::getline( infile, line, '\n' )) {
        std::istringstream line_ss(line);
        std::getline( line_ss, field1, ',' );
        std::getline( line_ss, field2, ',' );
        std::getline( line_ss, field3, ',' );
        std::getline( line_ss, field4, ',' );
        std::getline( line_ss, field5, ',' );
        std::getline( line_ss, field6, ',' );
        int j1,j2;
        long k1,k2;

        j1 = atoi(field2.c_str());
        k1 = atol(field3.c_str());
        j2 = atoi(field5.c_str());
        k2 = atol(field6.c_str());

        if (strcmp(field1.c_str(),"wavelet")==0 && strcmp(field4.c_str(),"wavelet")==0) {
            Index1D index_x(j1,k1,XWavelet);
            Index1D index_y(j2,k2,XWavelet);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else if (strcmp(field1.c_str(),"wavelet")==0 && strcmp(field4.c_str(),"scaling")==0) {
            Index1D index_x(j1,k1,XWavelet);
            Index1D index_y(j2,k2,XBSpline);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else if (strcmp(field1.c_str(),"scaling")==0 && strcmp(field4.c_str(),"wavelet")==0) {
            Index1D index_x(j1,k1,XBSpline);
            Index1D index_y(j2,k2,XWavelet);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else if (strcmp(field1.c_str(),"scaling")==0 && strcmp(field4.c_str(),"scaling")==0) {
            Index1D index_x(j1,k1,XBSpline);
            Index1D index_y(j2,k2,XBSpline);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else {
            std::cerr << "Got " << field1 << ", could not read file." << std::endl;
            exit(1); return;
        }
    }
}

void
getSparseGridIndexSet(const PrimalBasis &basis, IndexSet<Index2D> &Lambda, int j, T gamma)
{
    int j0 = basis.j0;
    for (long k1=basis.mra.rangeI(j0).firstIndex(); k1<=basis.mra.rangeI(j0).lastIndex(); ++k1) {
        for (long k2=basis.mra.rangeI(j0).firstIndex(); k2<=basis.mra.rangeI(j0).lastIndex(); ++k2) {
            Index1D row(j0,k1,XBSpline);
            Index1D col(j0,k2,XBSpline);
            Lambda.insert(Index2D(row,col));
        }
        for (int i2=1; i2<=j; ++i2) {
            int j2=j0+i2-1;
            for (long k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
                Index1D row(j0,k1,XBSpline);
                Index1D col(j2,k2,XWavelet);
                Lambda.insert(Index2D(row,col));
                Lambda.insert(Index2D(col,row));
            }
        }
    }
    for (int i1=1; i1<=j; ++i1) {
        int j1=j0+i1-1;
        for (int i2=1; i2<=j; ++i2) {
            if (T(i1+i2)-gamma*max(i1,i2)>(1-gamma)*j) {
                continue;
            }
            int j2=j0+i2-1;
            for (long k1=basis.rangeJ(j1).firstIndex(); k1<=basis.rangeJ(j1).lastIndex(); ++k1) {
                for (long k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
                    Index1D row(j1,k1,XWavelet);
                    Index1D col(j2,k2,XWavelet);
                    Lambda.insert(Index2D(row,col));
                }
            }
        }
    }
    return;
}

void
extendRHSIndexSet(const PrimalBasis &basis, IndexSet<Index2D> &Lambda, int J)
{
    std::cerr << "   extendRHSIndexSet adds level up to " << J << std::endl;
    IndexSet<Index1D> LambdaBoundary;
    for (int k=basis.mra.rangeI(basis.j0).firstIndex(); k<=basis.mra.rangeI(basis.j0).lastIndex(); ++k) {
        LambdaBoundary.insert(Index1D(basis.j0,k,XBSpline));
    }
    for (int j=basis.j0; j<=J; ++j) {
        for (int k=basis.rangeJL(j).firstIndex(); k<=basis.rangeJL(j).lastIndex(); ++k) {
            LambdaBoundary.insert(Index1D(j,k,XWavelet));
        }
        for (int k=basis.rangeJR(j).firstIndex(); k<=basis.rangeJR(j).lastIndex(); ++k) {
            LambdaBoundary.insert(Index1D(j,k,XWavelet));
        }
    }
    for (const_set1d_it it_x=LambdaBoundary.begin(); it_x!=LambdaBoundary.end(); ++it_x) {
        for (const_set1d_it it_y=LambdaBoundary.begin(); it_y!=LambdaBoundary.end(); ++it_y) {
            Lambda.insert(Index2D(*it_x,*it_y));
        }
    }
}

void
mv(LocalLaplaceIdentityOp2D &localLaplaceIdentityOp2D,
   LocalIdentityLaplaceOp2D &localIdentityLaplaceOp2D,
   Coefficients<Lexicographical,T,Index2D> &P,
   const IndexSet<Index2D> &Lambda, Coefficients<Lexicographical,T,Index2D> &v,
   Coefficients<Lexicographical,T,Index2D> &intermediate,
   Coefficients<Lexicographical,T,Index2D> &LIIAv, Coefficients<Lexicographical,T,Index2D> &IAUIv,
   Coefficients<Lexicographical,T,Index2D> &Av, T &time1, T &time2)
{
    //cerr << "      #v = " << v.size() << ", #Av = " << Lambda.size() << endl;
    Timer time;
    intermediate.setToZero();
    LIIAv.setToZero();
    FillWithZeros(Lambda,LIIAv);
    IAUIv.setToZero();
    FillWithZeros(Lambda,IAUIv);
    Av.setToZero();

    if (LIIAv.size()!=Lambda.size()) {
        LIIAv.clear();
        FillWithZeros(Lambda,LIIAv);
    }
    if (IAUIv.size()!=Lambda.size()) {
        IAUIv.clear();
        FillWithZeros(Lambda,IAUIv);
    }

    for (coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= P[(*it).first];
    }
    time.start();
    localLaplaceIdentityOp2D.evalAA(v,intermediate,LIIAv,IAUIv);
    time.stop();
    time1 = time.elapsed();
    //cerr << "      #LIIAv = " << LIIAv.size() << ", #IAUIv = " << IAUIv.size() << endl;

    Av += LIIAv;
    Av += IAUIv;
    LIIAv.setToZero();
    IAUIv.setToZero();
    time.start();
    localIdentityLaplaceOp2D.evalAA(v,intermediate,LIIAv,IAUIv);
    time.stop();
    time2 = time.elapsed();
    Av += LIIAv;
    Av += IAUIv;
    //cerr << "      #LIIAv = " << LIIAv.size() << ", #IAUIv = " << IAUIv.size() << endl;
    LIIAv.setToZero();
    IAUIv.setToZero();
    for (coeff2d_it it=Av.begin(); it!=Av.end(); ++it) {
        (*it).second *= P[(*it).first];
    }
    for (coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= 1./P[(*it).first];
    }
    std::cerr << "      MV: dof = " << Av.size() << ", time1 = " << time1
              << ", time2 = " << time2 << std::endl;
}

/*
T P_Lambda_r_norm_square = 0.;
for (const_set2d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
    P_Lambda_r_norm_square += std::pow(r[(*it)],2.);
    r.erase((*it));
}
T bulkalpha = 0.6;
T threshbound = std::sqrt(1-bulkalpha*bulkalpha) * r.norm(2.)/std::sqrt(T(r.size()));
Coefficients<Bucket,T,Index2D> r_bucket;
std::cerr << "      norm of r = " << r_norm << std::endl;
std::cerr << "      size of r = " << r.size() << std::endl;

r_bucket.bucketsort(r, threshbound);
IndexSet<Index2D> refinements;
for (int i=0; i<(int)r_bucket.bucket_ell2norms.size(); ++i) {
    P_Lambda_r_norm_square += std::pow(r_bucket.bucket_ell2norms[i],2.0L);
    int addDoF = r_bucket.addBucketToIndexSet(refinements,i);
    std::cerr << "      Added " << addDoF << " indices, now ||P_{Lambda}r ||_2 = "
              << std::sqrt(P_Lambda_r_norm_square)
              << ", alpha*r_norm = " << bulkalpha*r_norm << std::endl;
    if (P_Lambda_r_norm_square >= bulkalpha*r_norm*bulkalpha*r_norm) {
        int addDoF = r_bucket.addBucketToIndexSet(refinements,i+1);
        std::cerr << "      Added " << addDoF << " indices, now ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_r_norm_square) << std::endl;
        break;
    }
}

cout << " Refinements: " << refinements << endl;
for (const_set2d_it it=refinements.begin(); it!=refinements.end(); ++it) {
    extendMultiTree2(basis2d,(*it),offset,Lambda);
}
*/
/*
Coefficients<Lexicographical,T,Index1D> u_multi;
for (const_set2d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
    XType xtype1 = (*it).index1.xtype;
    int j1 = (*it).index1.j;
    long k1 = (*it).index1.k;
    if (xtype1==XBSpline && k1==basis.mra.rangeI(basis.j0).firstIndex()) {
        u_multi[(*it).index2] = 1.;
    }
}
stringstream filename2;
filename2 << "u_multi_" << example << "_" << d << "_" << threshTol << "_" << ell << "_" << iter << ".dat";
plotCoeff<T,PrimalBasis>(u_multi, basis, filename2.str().c_str(), false, true);
*/

/*
if (example==2) {
            for (coeff2d_it it=f.begin(); it!=f.end(); ++it) {
                XType xtype1 = (*it).first.index1.xtype;
                XType xtype2 = (*it).first.index2.xtype;
                int  j1 = (*it).first.index2.j;
                int  j2 = (*it).first.index2.j;
                long k1 = (*it).first.index2.k;
                long k2 = (*it).first.index2.k;
                if (xtype1==XWavelet && basis.rangeJI(j1).firstIndex()<=k1
                                     && k1<=basis.rangeJI(j1).lastIndex()) {
                    (*it).second=0.;
                }
                if (xtype2==XWavelet && basis.rangeJI(j2).firstIndex()<=k2
                                     && k2<=basis.rangeJI(j2).lastIndex()) {
                    (*it).second=0.;
                }
            }
        }

 */
