#include <iostream>
#include <lawa/lawa.h>
#include <boost/timer.hpp>

using namespace std;
using namespace lawa;

typedef double T;

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


typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef IndexSet<Index3D>::const_iterator                           const_set3d_it;
typedef Coefficients<Lexicographical,T,Index3D>::iterator           coeff3d_it;
typedef Coefficients<Lexicographical,T,Index3D>::const_iterator     const_coeff3d_it;

int main (int argc, char *argv[]) {

    cout.precision(20);
    if (argc!=4) {
        cout << "Usage: " << argv[0] << " d j0 J" << endl;
        return 0;
    }
    int example = 2;
    int d   = atoi(argv[1]);
    int j0  = atoi(argv[2]);
    int J  = atoi(argv[3]);
    T alpha = 0.7;
    T gamma = 0.005;
    const char* residualType = "standard";


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

    Timer time;
    stringstream performfilename;
    performfilename << "perform_multitree_mw_awgm_poisson3d_" << example << "_" << argv[1] << "_"
                  << argv[2] << "_" << alpha << "_" << gamma << "_" << residualType << ".dat";
    ofstream performfile(performfilename.str().c_str());
    Coefficients<Lexicographical,T,Index3D> u(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index3D> Au(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index3D> res(SIZEHASHINDEX2D);

    for (int iter=1; iter<=100; ++iter) {
        u.setToZero();
        Au.setToZero();
        res.setToZero();

        T time_galerkin1 = 0., time_galerkin2 = 0.;
        T time_residual1 = 0., time_residual2 = 0.;

        stringstream coefffilename;
        coefffilename << "coeff3d/coeff_multitree_mw_awgm_poisson3d_" << example << "_" << argv[1] << "_"
                      << argv[2] << "_" << alpha << "_" << gamma << "_" << residualType
                      << "__" << iter << ".dat";
        readCoefficientsFromFile(u, coefffilename.str().c_str());
        if (u.size()==0) break;

        Au = u;
        Au.setToZero();

        cout << "   Computing first matrix vector product for Galerkin system..." << endl;
        time.start();
        boost::timer boost_timer_galerkin1;
        localOp3D.eval(u,Au,Prec,"galerkin");
        double boost_time_galerkin1 = boost_timer_galerkin1.elapsed();
        time.stop();
        time_galerkin1 = time.elapsed();
        cout << "   ... finished." << endl;

        Au.setToZero();
        cout << "   Computing second matrix vector product for Galerkin system..." << endl;
        time.start();
        boost::timer boost_timer_galerkin2;
        localOp3D.eval(u,Au,Prec,"galerkin");
        double boost_time_galerkin2 = boost_timer_galerkin2.elapsed();
        time.stop();
        time_galerkin2 = time.elapsed();
        cout << "   ... finished." << endl;

        cout << "   Computing multi-tree for residual..." << endl;
        extendMultiTree(basis3d, u, res, residualType);
        cout << "   ... finished." << endl;

        cout << "   Computing first matrix vector product for residual..." << endl;
        time.start();
        boost::timer boost_timer_residual1;
        localOp3D.eval(u,res,Prec);
        double boost_time_residual1 = boost_timer_residual1.elapsed();
        time.stop();
        time_residual1 = time.elapsed();
        cout << "   ... finished." << endl;

        Au.setToZero();
        cout << "   Computing second matrix vector product for residual..." << endl;
        time.start();
        boost::timer boost_timer_residual2;
        localOp3D.eval(u,res,Prec);
        double boost_time_residual2 = boost_timer_residual2.elapsed();
        time.stop();
        time_residual2 = time.elapsed();
        cout << "   ... finished." << endl;

        cout << u.size() << " " << res.size() << " "
                         << time_galerkin1 << " " << time_galerkin2
                         << " " << time_residual1 << " " << time_residual2 << endl;
        cout << u.size() << " " << res.size() << " "
                         << boost_time_galerkin1 << " " << boost_time_galerkin2
                         << " " << boost_time_residual1 << " " << boost_time_residual2 << endl << endl;


        performfile << u.size() << " " << res.size()
                    << " " << time_galerkin1 << " " << time_galerkin2
                    << " " << time_residual1 << " " << time_residual2 << endl;
    }

    return 0;
}
