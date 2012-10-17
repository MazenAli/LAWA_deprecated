#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

typedef Basis<T,Orthogonal,Interval,Multi>                          PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;
typedef TensorBasis3D<Adaptive,PrimalBasis,PrimalBasis,PrimalBasis> Basis3D;

typedef OptimizedH1Preconditioner3D<T,Basis3D>                      Preconditioner3D;

///  Underlying bilinear form
typedef RefinementBasis::LaplaceOperator1D                          RefinementLaplaceOp1D;
typedef AdaptiveLaplaceOperator1D<T,Orthogonal,Interval,Multi>      LaplaceOp1D;

///  Local operator in 1d
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementLaplaceOp1D, LaplaceOp1D>         LocalOp1D;

typedef UniDirectionalLocalOperator<Index3D,XOne,LocalOp1D,
                                            NotXOne,Index2D>        UniDirectionalLocalOpXOne3D;
typedef UniDirectionalLocalOperator<Index3D,XTwo,LocalOp1D,
                                            NotXTwo,Index2D>        UniDirectionalLocalOpXTwo3D;
typedef UniDirectionalLocalOperator<Index3D,XThree,LocalOp1D,
                                            NotXThree,Index2D>      UniDirectionalLocalOpXThree3D;

typedef CompoundLocalOperator<Index3D,
                              UniDirectionalLocalOpXOne3D,
                              UniDirectionalLocalOpXTwo3D,
                              UniDirectionalLocalOpXThree3D>        CompoundLocalOperator3D;

//Righthandsides definitions (separable)
typedef RHSWithPeaks1D<T,PrimalBasis>                               Rhs1D;
typedef AdaptiveSeparableRhs<T,Index3D,Rhs1D,Rhs1D,Rhs1D>           AdaptiveSeparableRhsIntegral3D;
typedef CompoundRhs<T,Index3D,
                    AdaptiveSeparableRhsIntegral3D,
                    AdaptiveSeparableRhsIntegral3D,
                    AdaptiveSeparableRhsIntegral3D>                 CompoundRhsIntegral3D;

typedef MultiTreeAWGM<Index3D,Basis3D,CompoundLocalOperator3D,
                      CompoundRhsIntegral3D,Preconditioner3D>       MultiTreeAWGM3D;

typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef IndexSet<Index3D>::const_iterator                           const_set3d_it;
typedef Coefficients<Lexicographical,T,Index3D>::iterator           coeff3d_it;
typedef Coefficients<Lexicographical,T,Index3D>::const_iterator     const_coeff3d_it;

int example = 2;
T u1(T x)   {    return 1.; }
T u2(T y)   {    return 1.; }
T u3(T z)   {    return 1.; }
T du1(T x)  {    return 0.; }
T du2(T y)  {    return 0.; }
T du3(T z)  {    return 0.; }
T ddu1(T x) {    return -33-(T)(1.L/3.L); }
T ddu2(T y) {    return -33-(T)(1.L/3.L); }
T ddu3(T z) {    return -33-(T)(1.L/3.L); }

long double EnergyErrorSquared = 14.20158453089639L*14.20158453089639L;

/*
int example = 3;
T u1(T x)   {    return x*x*(1-x)*(1-x); }
T u2(T y)   {    return y*y*(1-y)*(1-y); }
T u3(T z)   {    return z*z*(1-z)*(1-z); }
T du1(T x)  {    return 2*x*(1-x)*(1-x)-2*x*x*(1-x); }
T du2(T y)  {    return 2*y*(1-y)*(1-y)-2*y*y*(1-y); }
T du3(T z)  {    return 2*z*(1-z)*(1-z)-2*z*z*(1-z); }
T ddu1(T x) {    return 2*(1-x)*(1-x) - 8*x*(1-x) + 2*x*x; }
T ddu2(T y) {    return 2*(1-y)*(1-y) - 8*y*(1-y) + 2*y*y; }
T ddu3(T z) {    return 2*(1-z)*(1-z) - 8*z*(1-z) + 2*z*z; }

long double EnergyErrorSquared = 3.*(1.L/630.L * 1.L/630.L * 2.L/105.L);
*/
T f1(T x)   {   return -ddu1(x); }

T f2(T y)   {   return -ddu2(y); }

T f3(T z)   {   return -ddu3(z); }

T sol(T x, T y, T z) {   return u1(x) * u2(y) * u3(z); }

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
    const char* treeType = "sparsetree";//"gradedtree";
    bool IsMW = true;
    bool compute_f_minus_Au_error = false;
    bool writeCoefficientsToFile = false;
    T eps   = 1e-5;
    Timer time;

    /// Basis initialization
    //PrimalBasis       basis(d,d_,j0);
    PrimalBasis       basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis  &refinementbasis = basis.refinementbasis;
    Basis3D basis3d(basis,basis,basis);

    /// Operator initialization
    LaplaceOp1D                     laplaceOp1D(basis);
    LocalOp1D                       localOp1D(basis,basis,refinementbasis.LaplaceOp1D, laplaceOp1D);
    UniDirectionalLocalOpXOne3D     uniDirectionalOpXOne3D(localOp1D);
    //uniDirectionalOpXOne3D.setParameters(J, 49157, 6151);
    UniDirectionalLocalOpXTwo3D     uniDirectionalOpXTwo3D(localOp1D);
    //uniDirectionalOpXTwo3D.setParameters(J, 49157, 6151);
    UniDirectionalLocalOpXThree3D   uniDirectionalOpXThree3D(localOp1D);
    //uniDirectionalOpXThree3D.setParameters(J, 49157, 6151);
    CompoundLocalOperator3D         localOp3D(uniDirectionalOpXOne3D,uniDirectionalOpXTwo3D,
                                              uniDirectionalOpXThree3D);

    /// Initialization of preconditioner
    Preconditioner3D  Prec(basis3d,1.,1.,1.,0.);

    /// Initialization of rhs
    DenseVectorT sing_pts_x, sing_pts_y, sing_pts_z;
    DenseMatrixT no_deltas, deltas_x, deltas_y, deltas_z;
    int order = 20;
    if (example == 2) order = 2*d;
    if (example == 3) order = 4+2*d;

    Function<T>                    fct_u1(u1,sing_pts_x), fct_f1(f1,sing_pts_x);
    Function<T>                    fct_u2(u1,sing_pts_y), fct_f2(f2,sing_pts_y);
    Function<T>                    fct_u3(u3,sing_pts_z), fct_f3(f3,sing_pts_y);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u1(basis, fct_u1, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f1(basis, fct_f1, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u2(basis, fct_u2, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f2(basis, fct_f2, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u3(basis, fct_u3, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f3(basis, fct_f3, no_deltas, order);
    Coefficients<Lexicographical,T,Index1D> rhs_u1_data(SIZEHASHINDEX1D),
                                            rhs_f1_data(SIZEHASHINDEX1D),
                                            rhs_u2_data(SIZEHASHINDEX1D),
                                            rhs_f2_data(SIZEHASHINDEX1D),
                                            rhs_u3_data(SIZEHASHINDEX1D),
                                            rhs_f3_data(SIZEHASHINDEX1D);
    AdaptiveSeparableRhsIntegral3D rhs1(rhs_f1, rhs_f1_data, rhs_u2, rhs_u2_data,
                                        rhs_u3, rhs_u3_data);
    AdaptiveSeparableRhsIntegral3D rhs2(rhs_u1, rhs_u1_data, rhs_f2, rhs_f2_data,
                                        rhs_u3, rhs_u3_data);
    AdaptiveSeparableRhsIntegral3D rhs3(rhs_u1, rhs_u1_data, rhs_u2, rhs_u2_data,
                                        rhs_f3, rhs_f3_data);
    CompoundRhsIntegral3D          F(rhs1,rhs2,rhs3);

    /// Initialization of multi tree based adaptive wavelet Galerkin method
    MultiTreeAWGM3D multiTreeAWGM3D(basis3d, localOp3D, F, Prec);
    multiTreeAWGM3D.setParameters(alpha, gamma, residualType, treeType, IsMW, writeCoefficientsToFile);

    Coefficients<Lexicographical,T,Index3D> u(SIZEHASHINDEX2D);
    getSparseGridVector(basis3d,u,0,(T)0.2);

    stringstream convfilename;
    convfilename << "conv_multitree_mw_awgm_poisson3d_" << example << "_" <<argv[1] << "_"
                 << argv[2] << "_" << alpha << "_" << gamma << "_" << residualType << "_"
                 << treeType << ".dat";
    stringstream coefffilename;
    coefffilename << "coeff_multitree_mw_awgm_poisson3d_" << example << "_" << argv[1] << "_"
                 << argv[2] << "_" << alpha << "_" << gamma << "_" << residualType << "_" << treeType;
    multiTreeAWGM3D.cg_solve(u, eps, 100, 1e-2, EnergyErrorSquared,
                             convfilename.str().c_str(), coefffilename.str().c_str());
    /*
    for (int j=0; j<=20; ++j) {
        Coefficients<Lexicographical,T,Index3D> u(SIZEHASHINDEX2D);
        getSparseGridVector(basis3d,u,j,(T)0.2);
        multiTreeAWGM3D.cg_solve(u, eps, "convfile.txt", 1, EnergyErrorSquared);

    }
    */
    return 0;
}

/*
    Index1D index1d_1(8,1535,XWavelet);
    Index1D index1d_2(0,1,XBSpline);
    Index3D index1(index1d_1,index1d_1,index1d_2);
    Index3D index2(index1d_1,index1d_2,index1d_1);
    Index3D index3(index1d_2,index1d_1,index1d_1);
    Index1D index1d_3(5,191,XWavelet);
    Index1D index1d_4(6,383,XWavelet);
    Index3D index4(index1d_3,index1d_3,index1d_4);
    Index3D index5(index1d_3,index1d_4,index1d_3);
    Index3D index6(index1d_4,index1d_3,index1d_3);

    cout << "# supp u = " << u.size() << endl;
    completeMultiTree(basis3d, index1, u);
    cout << "# supp u = " << u.size() << endl;
    completeMultiTree(basis3d, index2, u);
    cout << "# supp u = " << u.size() << endl;
    completeMultiTree(basis3d, index3, u);
    cout << "# supp u = " << u.size() << endl;
    completeMultiTree(basis3d, index4, u);
    cout << "# supp u = " << u.size() << endl;
    completeMultiTree(basis3d, index5, u);
    cout << "# supp u = " << u.size() << endl;
    completeMultiTree(basis3d, index6, u);
    cout << "# supp u = " << u.size() << endl;

*/

/*
 * IntegralF<Gauss,PrimalBasis> integral_u1(fct_u1,basis); integral_u1.quadrature.setOrder(order);
    IntegralF<Gauss,PrimalBasis> integral_u2(fct_u2,basis); integral_u2.quadrature.setOrder(order);
    IntegralF<Gauss,PrimalBasis> integral_u3(fct_u3,basis); integral_u3.quadrature.setOrder(order);
    IntegralF<Gauss,PrimalBasis> integral_f1(fct_f1,basis); integral_f1.quadrature.setOrder(order);
    IntegralF<Gauss,PrimalBasis> integral_f2(fct_f2,basis); integral_f2.quadrature.setOrder(order);
    IntegralF<Gauss,PrimalBasis> integral_f3(fct_f3,basis); integral_f3.quadrature.setOrder(order);
    for (int j1=j0+1; j1<=j0+1; ++j1) {
    for (int k1=basis.rangeJ(j1).firstIndex(); k1<=basis.rangeJ(j1).lastIndex(); ++k1) {
        Index1D index_x(j1,k1,XWavelet);
        for (int j2=j0+1; j2<=j0+1; ++j2) {
        for (int k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
            Index1D index_y(j2,k2,XWavelet);
            for (int j3=j0+1; j3<=j0+1; ++j3) {
            for (int k3=basis.rangeJ(j3).firstIndex(); k3<=basis.rangeJ(j3).lastIndex(); ++k3) {
                Index1D index_z(j3,k3,XWavelet);
                Index3D index(index_x,index_y,index_z);
                T val1 = integral_f1(j1,k1,XWavelet,0) * integral_u2(j2,k2,XWavelet,0) * integral_u3(j3,k3,XWavelet,0);
                T val2 = rhs1(index);
                T val3 = integral_u1(j1,k1,XWavelet,0) * integral_f2(j2,k2,XWavelet,0) * integral_u3(j3,k3,XWavelet,0);
                T val4 = rhs2(index);
                T val5 = integral_u1(j1,k1,XWavelet,0) * integral_u2(j2,k2,XWavelet,0) * integral_f3(j3,k3,XWavelet,0);
                T val6 = rhs3(index);
                if (fabs(val1-val2)>1e-16) cout << "Error integral 1 for " << index << endl;
                if (fabs(val3-val4)>1e-16) cout << "Error integral 2 for " << index << endl;
                if (fabs(val5-val6)>1e-16) cout << "Error integral 3 for " << index << endl;
                T val7 = val1 + val3 + val5;
                T val8 = F(index);
                if (fabs(val7-val8)>1e-16) cout << "Error integral 4 for " << index << endl;
            }
            }
        }
        }
    }
    }
    return 0;
 */
