#include <iostream>
#include <lawa/lawa.h>
#include <applications/canutoproject/lineartensorinterpolationpic2d.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

typedef Basis<T,Orthogonal,Interval,Multi>                          PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;
typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>             Basis2D;

typedef OptimizedH1Preconditioner2D<T,Basis2D>                      Preconditioner2D;
typedef NoPreconditioner<T,Index2D>                                 NoPreconditioner2D;

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
                              UniDirectionalLocalOpXTwo2D>           CompoundLocalOperator2D;
typedef ThetaTimeStepLocalOperator<Index2D, CompoundLocalOperator2D> ThetaTimeStepLocalOperator2D;

//Righthandsides definitions (separable)
typedef SmoothRHSWithAlignedSing2D<T, Basis2D, FullGridGL>          Integral_LinearTensorInterpolPic2D;
typedef RHS2D<T, Integral_LinearTensorInterpolPic2D,
              NoPreconditioner2D >                                  Rhs_LinearTensorInterpolPic2D;
typedef CompoundRhs<T,Index2D,Rhs_LinearTensorInterpolPic2D,
                    Rhs_LinearTensorInterpolPic2D,
                    Rhs_LinearTensorInterpolPic2D>                  CompoundRhs_LinearTensorInterpolPic2D;

typedef MultiTreeAWGM<Index2D,Basis2D,ThetaTimeStepLocalOperator2D,
                      CompoundRhs_LinearTensorInterpolPic2D,
                      Preconditioner2D>                             MultiTreeAWGM2D;

typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::iterator           coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;



int main (int argc, char *argv[]) {

    cout.precision(20);

    int d   = 2;
    int j0  = 0;
    T alpha = 0.7;
    T gamma = 0.001;
    const char* residualType = "standard";
    const char* treeType = "sparsetree"; //"gradedtree";
    bool IsMW = true;
    bool compute_f_minus_Au_error = false;
    bool writeCoefficientsToFile = true;
    T eps   = 1e-5;
    Timer time;


    /// Basis initialization
    //PrimalBasis       basis(d,d_,j0);
    PrimalBasis       basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis  &refinementbasis = basis.refinementbasis;
    Basis2D basis2d(basis,basis);

    /// Operator initialization
    LaplaceOp1D                  laplaceOp1D(basis);
    LocalOp1D                    localOp1D(basis,basis,refinementbasis.LaplaceOp1D,laplaceOp1D);
    UniDirectionalLocalOpXOne2D  uniDirectionalOpXOne2D(localOp1D,1.);
    UniDirectionalLocalOpXTwo2D  uniDirectionalOpXTwo2D(localOp1D,1.);
    CompoundLocalOperator2D      localOp2D(uniDirectionalOpXOne2D,uniDirectionalOpXTwo2D);
    ThetaTimeStepLocalOperator2D localOpPlusId2D(1.,1.,localOp2D);

    /// Initialization of preconditioner
    Preconditioner2D  Prec(basis2d,1.,1.,1.);
    NoPreconditioner2D  NoPrec;

    /// Initialization of rhs
    LinearTensorInterpolationPic2D<T> linearTensorInterpolPic2D;
    cout << "Reading image data and computing linear tensor interpolation of image data..." << endl;
    linearTensorInterpolPic2D.readPicture("claudiocanutogrey2.txt");
    cout << "... finished." << endl;
    Function2D<T> func_dx1_linearTensorInterpolPic2D(LinearTensorInterpolationPic2D<T>::dx1_evaluateInterpolationMinusLiftingFunction,
                                                     LinearTensorInterpolationPic2D<T>::sing_pts_x,
                                                     LinearTensorInterpolationPic2D<T>::sing_pts_y);
    Function2D<T> func_dx2_linearTensorInterpolPic2D(LinearTensorInterpolationPic2D<T>::dx2_evaluateInterpolationMinusLiftingFunction,
                                                     LinearTensorInterpolationPic2D<T>::sing_pts_x,
                                                     LinearTensorInterpolationPic2D<T>::sing_pts_y);
    Function2D<T> func_linearTensorInterpolPic2D(LinearTensorInterpolationPic2D<T>::evaluateInterpolationMinusLiftingFunction,
                                                 LinearTensorInterpolationPic2D<T>::sing_pts_x,
                                                 LinearTensorInterpolationPic2D<T>::sing_pts_y);
    int order = 1;
    Integral_LinearTensorInterpolPic2D      integral1_linearTensorInterpolPic2D(basis2d, func_dx1_linearTensorInterpolPic2D, order, 1, 0);
    Integral_LinearTensorInterpolPic2D      integral2_linearTensorInterpolPic2D(basis2d, func_dx2_linearTensorInterpolPic2D, order, 0, 1);
    Integral_LinearTensorInterpolPic2D      integral3_linearTensorInterpolPic2D(basis2d, func_linearTensorInterpolPic2D, order, 0, 0);
    Rhs_LinearTensorInterpolPic2D           rhs1_linearTensorInterpolPic2D(integral1_linearTensorInterpolPic2D, NoPrec);
    Rhs_LinearTensorInterpolPic2D           rhs2_linearTensorInterpolPic2D(integral2_linearTensorInterpolPic2D, NoPrec);
    Rhs_LinearTensorInterpolPic2D           rhs3_linearTensorInterpolPic2D(integral3_linearTensorInterpolPic2D, NoPrec);
    CompoundRhs_LinearTensorInterpolPic2D   F_linearTensorInterpolPic2D(rhs1_linearTensorInterpolPic2D,
                                                                        rhs2_linearTensorInterpolPic2D,
                                                                        rhs3_linearTensorInterpolPic2D);

    /// Initialization of multi tree based adaptive wavelet Galerkin method
    MultiTreeAWGM2D multiTreeAWGM2D(basis2d, localOpPlusId2D, F_linearTensorInterpolPic2D, Prec);
    multiTreeAWGM2D.setParameters(alpha, gamma, residualType, treeType, IsMW,
                                  writeCoefficientsToFile);

    Coefficients<Lexicographical,T,Index2D> u(SIZEHASHINDEX2D), f(SIZEHASHINDEX2D), Au(SIZEHASHINDEX2D);
    for (int k1=basis.mra.rangeI(j0).firstIndex(); k1<=basis.mra.rangeI(j0).lastIndex(); ++k1) {
        Index1D index1(j0,k1,XBSpline);
        for (int k2=basis.mra.rangeI(j0).firstIndex(); k2<=basis.mra.rangeI(j0).lastIndex(); ++k2) {
            Index1D index2(j0,k2,XBSpline);
            Index2D index(index1,index2);
            u[index] = integral3_linearTensorInterpolPic2D(index);
            f[index] = Prec(index) * F_linearTensorInterpolPic2D(index);
            Au[index] = 0.;
        }
    }

    localOpPlusId2D.eval(u, Au, Prec, "galerkin");

    cout << "u = " << u << endl;
    cout << "Au = " << Au << endl;
    cout << "f = " <<  f << endl;

    Au -= f;
    cout << "Error norm: " << Au.norm(2.) << endl;


    stringstream convfilename;
    convfilename << "conv_poisson2d_pic_" << alpha << "_" << gamma << "_" << residualType << "_"
                 << treeType << ".dat";
    stringstream coefffilename;
    coefffilename << "coeff_poisson2d_pic_" << alpha << "_" << gamma << "_" << residualType << "_"
                  << treeType;

    multiTreeAWGM2D.cg_solve(u, eps, 1, 1e-2, 0., convfilename.str().c_str(),
                             coefffilename.str().c_str());

//    plot2D<T,Basis2D,Preconditioner>(basis2d, u, Prec, sol, 0., 1., 0., 1., 0.1, "multiTreeAWGM_sol");

    return 0;
}
