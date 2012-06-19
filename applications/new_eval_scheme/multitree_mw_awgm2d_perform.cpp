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

typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::iterator           coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

//void
//readCoefficientsFromFile(const PrimalBasis &basis, Coefficients<Lexicographical,T,Index2D> &u, const char* filename);

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
    T eps   = 1e-9;

    /// Basis initialization
    PrimalBasis       basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis  &refinementbasis = basis.refinementbasis;
    Basis2D basis2d(basis,basis);

    /// Operator initialization
    LaplaceOp1D                  laplaceOp1D(basis);
    LocalOp1D                    localOp1D(basis,basis,refinementbasis.LaplaceOp1D,laplaceOp1D);
    UniDirectionalLocalOpXOne2D  uniDirectionalOpXOne2D(localOp1D);
    uniDirectionalOpXOne2D.setParameters(8, 6151, 193);
    UniDirectionalLocalOpXTwo2D  uniDirectionalOpXTwo2D(localOp1D);
    uniDirectionalOpXOne2D.setParameters(8, 6151, 193);
    CompoundLocalOperator2D      localOp2D(uniDirectionalOpXOne2D,uniDirectionalOpXTwo2D);

    /// Initialization of preconditioner
    Preconditioner  Prec(basis2d,1.,1.,0.);


    Timer time;
    stringstream performfilename;
    performfilename << "perform_multitree_mw_awgm_poisson2d_" << example << "_" << argv[1] << "_"
                  << argv[2] << "_" << alpha << "_" << gamma << "_" << residualType << ".dat";
    ofstream performfile(performfilename.str().c_str());
    Coefficients<Lexicographical,T,Index2D> u(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> Au1(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> Au2(SIZEHASHINDEX2D);

    for (int iter=1; iter<=100; ++iter) {
        u.setToZero();
        Au1.setToZero();
        Au2.setToZero();
        T time_galerkin1 = 0., time_galerkin2 = 0.;
        T time_residual1 = 0., time_residual2 = 0.;
        stringstream coefffilename;
        coefffilename << "coeff2d/coeff_multitree_mw_awgm_poisson2d_" << example << "_" << argv[1] << "_"
                      << argv[2] << "_" << alpha << "_" << gamma << "_" << residualType
                      << "__" << iter << ".dat";
        readCoefficientsFromFile(u, coefffilename.str().c_str());
        if (u.size()==0) break;
/*
        Au1 = u;
        Au1.setToZero();

        time.start();
        localOp2D.eval(u,Au1,Prec);
        time.stop();
        time_galerkin1 = time.elapsed();
        Au1.setToZero();
        time.start();
        localOp2D.eval(u,Au1,Prec);
        time.stop();
        time_galerkin2 = time.elapsed();

        extendMultiTree(basis2d, u, Au2, residualType);
        time.start();
        localOp2D.eval(u,Au2,Prec);
        time.stop();
        time_residual1 = time.elapsed();
        Au2.setToZero();
        time.start();
        localOp2D.eval(u,Au2,Prec);
        time.stop();
        time_residual2 = time.elapsed();

        performfile << u.size() << " " << time_galerkin1 << " " << time_galerkin2
                                << " " << time_residual1 << " " << time_residual2 << endl;
*/
    }
    return 0;
}

/*
void
readCoefficientsFromFile(const PrimalBasis &basis, Coefficients<Lexicographical,T,Index2D> &u, const char* filename)
{
    std::ifstream infile (filename);
    if (infile.is_open()) {
        std::cout << "File is open, ready to read..." << std::endl;
        std::string line;
        Coefficients<Lexicographical,T,Index1D> tmp;
        while(std::getline( infile, line, '\n' )) {
            //cout << line << endl;
            std::string field1, field2, field3, field4, field5, field6, field7;
            std::istringstream line_ss(line);
            std::getline( line_ss, field1, ',' );
            std::getline( line_ss, field2, ',' );
            std::getline( line_ss, field3, ',' );
            std::getline( line_ss, field4, ',' );
            std::getline( line_ss, field5, ',' );
            std::getline( line_ss, field6, ' ' );
            std::getline( line_ss, field7, ' ' );
            //std::cerr << field1 << " " << field2 << " " << field3 << " " << field4
            //          << " " << field5 << " " << field6 << " " << field7 << std::endl;
            Index1D index1, index2;
            double val;
            int j1 = atoi(field2.c_str());
            int k1 = atoi(field3.c_str());
            index1.j = j1; index1.k = k1;
            int j2 = atoi(field5.c_str());
            int k2 = atoi(field6.c_str());
            index2.j = j2; index2.k = k2;
            val = atof(field7.c_str());
            //std::cerr << field1 << " " << j1 << " " << k1 << " " << field4 << " " << j2 << " " << k2 << " " << val << " " << field7 << endl << endl;

            if (strcmp(field1.c_str(),"wavelet")==0)    index1.xtype = XWavelet;
            else                                        index1.xtype = XBSpline;
            if (strcmp(field4.c_str(),"wavelet")==0)    index2.xtype = XWavelet;
            else                                        index2.xtype = XBSpline;

            Index2D index(index1,index2);
            u[index] = val;
            //if (index1.xtype==XBSpline && index1.k==1) {
            //    tmp[index2] = val;
            //}
        }
        //cout << "Size of 1d-tree: " << tmp.size() << endl;
        //plotCoeff(tmp, basis, "coeff", false, true);

    }
    else {
        std::cout << "File not found." << std::endl;
        exit(1);
        return;
    }
    //cout << "u = " << u << endl;
}
*/
