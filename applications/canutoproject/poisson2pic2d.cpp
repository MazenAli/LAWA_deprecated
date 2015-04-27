#include <iostream>
#include <lawa/lawa.h>
#include <applications/canutoproject/lineartensorinterpolationpic2d.h>
#include <applications/canutoproject/rhstensorinterpolationpic2d.h>

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
typedef RHSTensorInterpolationPic2D<Basis2D,
                                    ThetaTimeStepLocalOperator2D>    RhsTensorInterpolationPic2D;

typedef MultiTreeAWGM<Index2D,Basis2D,ThetaTimeStepLocalOperator2D,
                      RhsTensorInterpolationPic2D,
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
    int j0  = 3;
    T alpha = 0.7;
    T gamma = 0.001;
    const char* residualType = "standard";
    const char* treeType = "sparsetree"; //"gradedtree";
    bool IsMW = true;
    bool compute_f_minus_Au_error = false;
    bool writeCoefficientsToFile = true;
    T eps   = 1e-5;
    Timer time;


    T diffusion_coeff = 0.000001;

    /// Basis initialization
    //PrimalBasis       basis(d,d_,j0);
    PrimalBasis       basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis  &refinementbasis = basis.refinementbasis;
    Basis2D basis2d(basis,basis);

    /// Operator initialization
    LaplaceOp1D                  laplaceOp1D(basis);
    LocalOp1D                    localOp1D(basis,basis,refinementbasis.LaplaceOp1D,laplaceOp1D);
    UniDirectionalLocalOpXOne2D  uniDirectionalOpXOne2D(localOp1D,diffusion_coeff);
    UniDirectionalLocalOpXTwo2D  uniDirectionalOpXTwo2D(localOp1D,diffusion_coeff);
    CompoundLocalOperator2D      localOp2D(uniDirectionalOpXOne2D,uniDirectionalOpXTwo2D);
    ThetaTimeStepLocalOperator2D localOpPlusId2D(1.,1.,localOp2D);

    /// Initialization of preconditioner
    Preconditioner2D  Prec(basis2d,diffusion_coeff,diffusion_coeff,1.);
    NoPreconditioner2D  NoPrec;

    /// Initialization of rhs
    LinearTensorInterpolationPic2D<T> linearTensorInterpolPic2D;
    cout << "Reading image data and computing linear tensor interpolation of image data..." << endl;
    linearTensorInterpolPic2D.readPicture("claudiocanutogrey3.txt");
    cout << "... finished." << endl;

    Coefficients<Lexicographical,T,Index2D> u_ref(SIZEHASHINDEX2D), f_init;
    //readCoefficientsFromFile(u_ref,"claudiocanutogrey2_L2/coeff_image_L2_0.7_0.1_standard_sparsetree__15.dat");
    readCoefficientsFromFile(u_ref,"claudiocanutogrey3_L2/coeff_image_L2_0.7_0.1_standard_sparsetree__16.dat");
    cout << "Size of reference solution vector: #supp u_ref = " << u_ref.size() << endl;
    getSparseGridVector(basis2d, f_init, 7, 0.);
    f_init += u_ref;
    f_init.setToZero();
    cout << "Size of initial rhs vector: #supp f_init = " << f_init.size() << endl;
    RhsTensorInterpolationPic2D rhsTensorInterpolationPic2D(basis2d, localOpPlusId2D, u_ref);
    rhsTensorInterpolationPic2D.initializeRHS(f_init);

    /// Initialization of multi tree based adaptive wavelet Galerkin method
    MultiTreeAWGM2D multiTreeAWGM2D(basis2d, localOpPlusId2D, rhsTensorInterpolationPic2D, Prec);
    multiTreeAWGM2D.setParameters(alpha, gamma, residualType, treeType, IsMW,
                                  writeCoefficientsToFile);

    Coefficients<Lexicographical,T,Index2D> u(SIZEHASHINDEX2D);
    for (int k1=basis.mra.rangeI(j0).firstIndex(); k1<=basis.mra.rangeI(j0).lastIndex(); ++k1) {
        Index1D index1(j0,k1,XBSpline);
        for (int k2=basis.mra.rangeI(j0).firstIndex(); k2<=basis.mra.rangeI(j0).lastIndex(); ++k2) {
            Index1D index2(j0,k2,XBSpline);
            Index2D index(index1,index2);
            u[index] = 0.;
        }
    }

    stringstream convfilename;
    convfilename << "conv_poisson2pic2d_" << alpha << "_" << gamma << "_" << diffusion_coeff
                 << "_" << residualType << "_" << treeType << ".dat";
    stringstream coefffilename;
    coefffilename << "coeff_poissonpic2pic2d_" << alpha << "_" << gamma << "_" << diffusion_coeff
                  << "_" << residualType << "_" << treeType;

    //multiTreeAWGM2D.cg_solve(u, eps, 100, 1e-2, 0., convfilename.str().c_str(),
    //                         coefffilename.str().c_str());
    cerr << "Warning: cg solver not started, RHS type incompatible with cg_solve, no propagation present" << endl;

    return 0;
}

