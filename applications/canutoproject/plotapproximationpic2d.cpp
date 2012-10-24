#include <iostream>
#include <fstream>
#include <lawa/lawa.h>
#include <applications/canutoproject/lineartensorinterpolationpic2d.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;

typedef Basis<T,Orthogonal,Interval,Multi>                          PrimalBasis;
typedef TensorBasis2D<Adaptive, PrimalBasis, PrimalBasis>           Basis2D;

typedef OptimizedH1Preconditioner2D<T,Basis2D>                      Preconditioner;

typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;

void
plotImageApproximation2D(const LinearTensorInterpolationPic2D<T> &linearTensorInterpolationPic,
                         const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D> &coeff,
                         DenseMatrixT &currentEvaluations, const char* filename,
                         T &L2Error, T &LinftyError);

void
plotImageScatterCoeff2D(const LinearTensorInterpolationPic2D<T> &linearTensorInterpolationPic,
                        const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D> &coeff,
                        const char* filename);



int main (int argc, char *argv[]) {

    int j0 = 3;
    PrimalBasis basis(2,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    Basis2D basis2d(basis, basis);

    cout.precision(6);
    if (argc!=3) {
        cout << "Usage: " << argv[0] << " approxType picturenumber" << endl;
        return 0;
    }
    int picturenumber = atoi(argv[2]);  // 2 = medium quality, 3 = high quality picture of C. Canuto

    T alpha  = 0.7;
    T gamma  = 0.001;
    const char* residualType = "standard";
    const char* treeType = "sparsetree"; //"gradedtree";

    T diffusion_coeff = 0.000001;
    Preconditioner  Prec(basis2d,diffusion_coeff,diffusion_coeff,1.);


    LinearTensorInterpolationPic2D<T> linearTensorInterpolPic2D;
    cout << "Reading image data and computing linear tensor interpolation of image data..." << endl;
    stringstream inputfilename;
    inputfilename << "claudiocanutogrey" << picturenumber << ".txt";
    linearTensorInterpolPic2D.readPicture(inputfilename.str().c_str());
    cout << "... finished." << endl;

    int PlotPts1 = linearTensorInterpolPic2D.N1+1, PlotPts2 = linearTensorInterpolPic2D.N2+1;
    cout << "PlotPts1 = " << PlotPts1 << ", PlotsPts2 = " << PlotPts2 << endl;

    size_t hashMapSize = 98317; // 196613
    Coefficients<Lexicographical,T,Index2D> u(hashMapSize), u_old(hashMapSize), u_leafs(hashMapSize);

    DenseMatrixT currentL2Evaluations(_(0,PlotPts1-1), _(0,PlotPts2-1));

    stringstream convfilename;
    convfilename << "conv_image_" << alpha << "_" << gamma << "_" << diffusion_coeff << "_" << residualType << "_" << treeType << ".dat";
    ofstream convfile(convfilename.str().c_str());
    for (int iter=1; iter<=100; ++iter) {
        T L2Error = 0., LinftyError = 0.;
        int N = 0;
        if (strcmp("L2Approximation",argv[1])==0) {
            stringstream coefffilename;
            coefffilename << "claudiocanutogrey" << picturenumber << "_L2/coeff_image_L2_"
                          << alpha << "_" << gamma << "_" << residualType << "_" << treeType
                          << "__" << iter << ".dat";
            readCoefficientsFromFile(u, coefffilename.str().c_str());
            N = u.size();
            u_leafs.clear();
            for (const_coeff2d_it it=u.begin(); it!=u.end(); ++it) {
                if (u_old.find((*it).first)==u_old.end()) u_leafs[(*it).first] = (*it).second;
            }

            stringstream plotfilename, scattercoefffilename;
            plotfilename << "image_L2_" << alpha << "_" << gamma << "_" << residualType << "_"
                         << treeType << "_" << iter;
            scattercoefffilename << "scattercoeff_image_L2_" << alpha << "_" << gamma << "_"
                                 << residualType << "_" << treeType << "_" << iter;

            plotImageScatterCoeff2D(linearTensorInterpolPic2D, basis2d, u, scattercoefffilename.str().c_str());

            L2Error = 0., LinftyError = 0.;
            if (iter==1) {
                plotImageApproximation2D(linearTensorInterpolPic2D, basis2d, u, currentL2Evaluations, plotfilename.str().c_str(), L2Error, LinftyError);
            }
            else {
                plotImageApproximation2D(linearTensorInterpolPic2D, basis2d, u_leafs, currentL2Evaluations, plotfilename.str().c_str(), L2Error, LinftyError);
            }

            u_old.clear();
            u_old = u;
        }
        else if (strcmp("PoissonApproximation",argv[1])==0) {
            // Coefficients are solutions of linear system and are __updated__ within each iteration
            // of the adaptive scheme. So, we cannot use former plot data.
            DenseMatrixT currentPoissonEvaluations(_(0,PlotPts1-1), _(0,PlotPts2-1));
            stringstream coefffilename;
            u.clear();
            coefffilename << "claudiocanutogrey" << picturenumber << "_poisson/coeff_poissonpic2pic2d_"
                          << alpha << "_" << gamma << "_" << diffusion_coeff << "_"
                          << residualType << "_" << treeType << "__" << iter << ".dat";
            readCoefficientsFromFile(u, coefffilename.str().c_str());
            N = u.size();

            for (coeff2d_it it=u.begin(); it!=u.end(); ++it) {
                (*it).second *= Prec((*it).first);
            }

            stringstream plotfilename, scattercoefffilename;
            plotfilename << "image_poisson_" << alpha << "_" << gamma << "_" << diffusion_coeff
                         << "_" << residualType << "_" << treeType << "_" << iter;
            scattercoefffilename << "scattercoeff_image_poisson_" << alpha << "_" << gamma << "_"
                                 << diffusion_coeff << "_" << residualType << "_" << treeType << "_" << iter;

            plotImageScatterCoeff2D(linearTensorInterpolPic2D, basis2d, u, scattercoefffilename.str().c_str());

            L2Error = 0., LinftyError = 0.;
            plotImageApproximation2D(linearTensorInterpolPic2D, basis2d, u, currentPoissonEvaluations, plotfilename.str().c_str(), L2Error, LinftyError);

        }
        else {
            std::cerr << "Unknown approximation type." << std::endl;
            exit(1);
        }
        convfile << iter << " " << N << " " << L2Error << " " << LinftyError << endl;
    }

    return 0.;
}


void
plotImageApproximation2D(const LinearTensorInterpolationPic2D<T> &linearTensorInterpolationPic,
                         const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D> &coeff,
                         DenseMatrixT &currentEvaluations,
                         const char* filename, T &L2Error, T &LinftyError)
{
    std::cerr << "Just for testing..." << std::endl;

    int PlotPts1 = currentEvaluations.numRows(), PlotPts2 = currentEvaluations.numCols();
    T h1 = 1./(PlotPts1-1), h2 = 1./(PlotPts2-1);
    std::stringstream PlotFileName;
    PlotFileName << filename << ".dat";
    std::ofstream plotfile(PlotFileName.str().c_str());
    plotfile.precision(8);

    for (const_coeff2d_it it = coeff.begin(); it != coeff.end(); ++it) {
        XType xtype1 = (*it).first.index1.xtype;
        XType xtype2 = (*it).first.index2.xtype;
        int j1 = (*it).first.index1.j, j2 = (*it).first.index2.j;
        long k1 = (*it).first.index1.k, k2 = (*it).first.index2.k;
        T coeff = (*it).second;

        Support<T> supp_x = basis.first.generator(xtype1).support(j1,k1);
        Support<T> supp_y = basis.second.generator(xtype2).support(j2,k2);

        int i1_first = floor(supp_x.l1 / h1), i1_last = ceil(supp_x.l2 /h1);
        int i2_first = floor(supp_y.l1 / h2), i2_last = ceil(supp_y.l2 /h2);

        for (T i1=i1_first; i1<=i1_last; ++i1) {
            T x1 = i1*h1;
            for (T i2=i2_first; i2<i2_last; ++i2) {
                T x2 = i2*h2;
                T tmp = coeff * basis.first.generator(xtype1)(x1,j1,k1,0) * basis.second.generator(xtype2)(x2,j2,k2,0);
                currentEvaluations(i1,i2) += tmp;
            }
        }

    }

    L2Error = 0.;
    LinftyError = 0.;
    for (int i1=0; i1<PlotPts1; ++i1) {
        T x1 = i1*h1;
        for (int i2=0; i2<PlotPts2; ++i2) {
            T x2 = i2*h2;

            T appr = currentEvaluations(i1,i2) + linearTensorInterpolationPic.evaluateLiftingFunction(x1,x2);
            T exact= linearTensorInterpolationPic.evaluateInterpolation(x1,x2);

            exact = std::max(exact,(T)0.); appr = std::max(appr,(T)0.);
            exact = std::min(exact,(T)1.); appr = std::min(appr,(T)1.);

            plotfile << i1 << " " << i2 << " " << exact << " " << appr  << std::endl;
            T localError = fabs(exact-appr);

            //if (localError == 1) {
            //    std::cerr << "(" << i1 << ", " << i2 << "): maximum error." << std::endl;
            //}

            LinftyError = std::max(LinftyError, fabs(localError));
            L2Error    += localError*localError;
        }
        plotfile << std::endl;
    }
    L2Error *= h1 * h2;
    L2Error = std::sqrt(L2Error);

    plotfile.close();

}

void
plotImageScatterCoeff2D(const LinearTensorInterpolationPic2D<T> &linearTensorInterpolationPic,
                        const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D> &coeff,
                        const char* filename)
{
    int N1 = linearTensorInterpolationPic.N1;
    int N2 = linearTensorInterpolationPic.N2;

    std::stringstream PlotFileName;
    PlotFileName << filename << ".dat";
    std::ofstream plotfile(PlotFileName.str().c_str());
    plotfile.precision(8);

    for (const_coeff2d_it it = coeff.begin(); it != coeff.end(); ++it) {
        int  j_x=(*it).first.index1.j, j_y=(*it).first.index2.j;
        long k_x=(*it).first.index1.k,  k_y=(*it).first.index2.k;
        XType xtype_x=(*it).first.index1.xtype, xtype_y=(*it).first.index2.xtype;

        Support<T> supp_x = basis.first.generator(xtype_x).support(j_x,k_x);
        Support<T> supp_y = basis.second.generator(xtype_y).support(j_y,k_y);

        T x=0., y=0.;
        x = N1*((supp_x.l2 + supp_x.l1)/(T)2.);
        y = N2*((supp_y.l2 + supp_y.l1)/(T)2.);

        plotfile << x << " " << y << " " << (*it).second << " " << -1. << std::endl;
    }
    plotfile.close();
    return;
}
