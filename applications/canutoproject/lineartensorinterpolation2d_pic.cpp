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

typedef Integral2D<FullGridGL, PrimalBasis, PrimalBasis>            Integral_LinearTensorInterpolPic2D_Test;
typedef SmoothRHSWithAlignedSing2D<T, Basis2D, FullGridGL>          Integral_LinearTensorInterpolPic2D;
typedef RHS2D<T, Integral_LinearTensorInterpolPic2D,
              NoPreconditioner<T,Index2D> >                         Rhs_LinearTensorInterpolPic2D;


typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;

//void
//plotImageApproximation2D(const LinearTensorInterpolationPic2D<T> linearTensorInterpolationPic,
//                         const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D> coeff,
//                         T a1, T b1, T a2, T b2, T h, const char* filename);

void
plotImageApproximation2D(const LinearTensorInterpolationPic2D<T> linearTensorInterpolationPic,
                         const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D> coeff,
                         DenseMatrixT &currentEvaluations, const char* filename, T &L2Error, T &LinftyError);

void
plotImageScatterCoeff2D(const LinearTensorInterpolationPic2D<T> linearTensorInterpolationPic,
                        const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D> coeff,
                        const char* filename);

int main (int argc, char *argv[]) {
    int j0 = 3;
    //PrimalBasis basis(2,2,0);
    PrimalBasis basis(2,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    Basis2D basis2d(basis, basis);

    LinearTensorInterpolationPic2D<T> linearTensorInterpolPic2D;
    cout << "Reading image data and computing linear tensor interpolation of image data..." << endl;
    linearTensorInterpolPic2D.readPicture("claudiocanutogrey3.txt");
    cout << "... finished." << endl;


    //cout << "Plotting interpolation of image data..." << endl;
    //linearTensorInterpolPic2D.plotInterpolation("claudiocanutogrey2_interpol.txt", 0.002, 0.002);
    //cout << "... finished." << endl;

    Function2D<T> func_linearTensorInterpolPic2D(LinearTensorInterpolationPic2D<T>::evaluateInterpolationMinusLiftingFunction,
                                                 LinearTensorInterpolationPic2D<T>::sing_pts_x,
                                                 LinearTensorInterpolationPic2D<T>::sing_pts_y);
    int order = 1, deriv_x = 0, deriv_y = 0;
    NoPreconditioner<T,Index2D> noPrec;

    Integral_LinearTensorInterpolPic2D integral_linearTensorInterpolPic2D(basis2d, func_linearTensorInterpolPic2D, order, deriv_x, deriv_y);
    Rhs_LinearTensorInterpolPic2D F(integral_linearTensorInterpolPic2D, noPrec);

    //int refinementLevel1 = 8, refinementLevel2 = 8;
    //int PlotPts1 = pow2i<T>(refinementLevel1)+1, PlotPts2 = pow2i<T>(refinementLevel2)+1;
    int PlotPts1 = linearTensorInterpolPic2D.N1+1, PlotPts2 = linearTensorInterpolPic2D.N2+1;
    cout << "PlotPts1 = " << PlotPts1 << ", PlotsPts2 = " << PlotPts2 << endl;
    DenseMatrixT currentEvaluations(_(0,PlotPts1-1), _(0,PlotPts2-1));

    //    ADAPTIVE SCHEME
    T alpha = 0.7;
    T gamma = 0.1;
    T eps   = 1e-8;

    const char* residualType = "standard";
    bool sparsetree = true;
    bool IsMW = true;
    int NumOfIterations = 100;

    size_t hashMapSize = 98317; // 196613
    Coefficients<Lexicographical,T,Index2D> u(hashMapSize), res(hashMapSize);
    Coefficients<Lexicographical,T,Index2D> p(hashMapSize);
    Coefficients<Lexicographical,T,Index2D> u_leafs(hashMapSize); // "leafs" of u


    for (int k1=basis.mra.rangeI(j0).firstIndex(); k1<=basis.mra.rangeI(j0).lastIndex(); ++k1) {
        Index1D index1(j0,k1,XBSpline);
        for (int k2=basis.mra.rangeI(j0).firstIndex(); k2<=basis.mra.rangeI(j0).lastIndex(); ++k2) {
            Index1D index2(j0,k2,XBSpline);
            Index2D index(index1,index2);
            u[index] = 0.;
        }
    }

    cout << "Size of initial u: " << u.size() << endl;
    for (const_coeff2d_it it=u.begin(); it!=u.end(); ++it) {
        u[(*it).first] = F((*it).first);
        u_leafs[(*it).first] = 0.;
        p[(*it).first] = 0.;
        res[(*it).first] = 0.;
    }

    stringstream convfilename;
    convfilename << "conv_image_" << alpha << "_" << gamma << "_" << residualType << ".dat";
    std::ofstream convfile(convfilename.str().c_str());

    T L2Error = 1., LinftyError = 1.;

    for (int iter=1; iter<=NumOfIterations; ++iter) {
        int N = u.size();
        int N_residual = 0;


        std::cerr << std::endl << "   *****   Iteration " << iter << " *****" << std::endl << std::endl;
        std::cerr << "      Current number of dof = " << u.size() << std::endl;

        // ------------------------------------------------------------------------
        // ----------------------- Resetting of vectors ---------------------------

        res.setToZero();
        p.clear();
        extendMultiTree(basis2d, u_leafs, res, residualType, IsMW, sparsetree);
        std::cerr << "      #supp u = " << u.size() << ", #supp r = " << res.size() << std::endl;
        std::cerr << "        Computing residual..." << std::endl;

        // Store current multi-tree structure of u
        u_leafs = u;

        for (coeff2d_it it=res.begin(); it!=res.end(); ++it) {
            if (u.find((*it).first)!=u.end()) {
                res.erase((*it).first);
            }
            else {
                (*it).second = F((*it).first);
            }
        }
        T Residual = res.norm(2.);

        convfile << u.size() << " " << Residual << " " << L2Error << " " << LinftyError << endl;

        std::cerr << "      ... finished with Residual: " << Residual << std::endl;
        if (Residual <= eps) {
            std::cerr << "      Target tolerance reached: Residual = " << Residual
                      << ", eps = " << eps << std::endl;
            break;
        }

        // ------------------------------------------------------------------------

        // --------------------- Computing next index set  ------------------------
        long double P_Lambda_Residual_square = 0.0L;
        if (res.size()!=0) {
            T threshbound = std::sqrt(1-alpha*alpha) * res.norm((T)2.)/std::sqrt(T(res.size()));
            Coefficients<Bucket,T,Index2D> r_bucket;
            r_bucket.bucketsort(res, threshbound);
            std::cerr << "         ||P_{Lambda}r ||_2 = " << std::sqrt(P_Lambda_Residual_square)
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
        for (const_coeff2d_it it=u.begin(); it!=u.end(); ++it) {
            res[(*it).first] = (*it).second;
        }

        // The vector p satisfies the bulk criterion. But it is not yet a multi-tree...
        std::cerr << "      Size of u before extension: " << u.size() << std::endl;
        for (const_coeff2d_it it=p.begin(); it!=p.end(); ++it) {
            if (u.find((*it).first)==u.end()) {
                completeMultiTree(basis2d, (*it).first, u, 0, sparsetree);
            }
        }
        std::cerr << "      Size of u after extension: " << u.size() << std::endl;

        // Note that u_leafs is still supported on the previous Galerkin index set!
        for (coeff2d_it it=u.begin(); it!=u.end(); ++it) {
            if (u_leafs.find((*it).first)!=u_leafs.end()) u_leafs.erase((*it).first);
            else {
                (*it).second = res[(*it).first];
                u_leafs[(*it).first] = res[(*it).first];;
            }
        }



        stringstream plotfilename, plotfilename2, coefffilename, scattercoefffilename;

        coefffilename  << "coeff_image_" << alpha << "_" << gamma << "_" << residualType;
        writeCoefficientsToFile(u,iter,coefffilename.str().c_str());

        plotfilename  << "image_test_" << alpha << "_" << gamma << "_" << residualType << "_" << iter;
        plotfilename2 << "image_" << alpha << "_" << gamma << "_" << residualType << "_" << iter;
        scattercoefffilename << "image_coeff_" << alpha << "_" << gamma << "_" << residualType << "_" << iter;

        plotImageScatterCoeff2D(linearTensorInterpolPic2D, basis2d, u, scattercoefffilename.str().c_str());

        //plotScatterCoeff(u, basis2d, coefffilename.str().c_str(), true);

        //plotImageApproximation2D(linearTensorInterpolPic2D, basis2d, u, 0., 1., 0., 1., pow2i<T>(-refinementLevel1), plotfilename.str().c_str());
        if (iter==1) {
            //plotImageApproximation2D(linearTensorInterpolPic2D, basis2d, u_leafs, currentEvaluations, 0., 1, 0., 1., plotfilename2.str().c_str());
            plotImageApproximation2D(linearTensorInterpolPic2D, basis2d, u, currentEvaluations, plotfilename2.str().c_str(), L2Error, LinftyError);
        }
        else {
            plotImageApproximation2D(linearTensorInterpolPic2D, basis2d, u_leafs, currentEvaluations, plotfilename2.str().c_str(), L2Error, LinftyError);
        }

        u_leafs.setToZero();

        // ------------------------------------------------------------------------
    }
    convfile.close();

    //plotfile.close();

    //plotLinearInterplolation(A, basis);

    return 0;

}


void
plotImageApproximation2D(const LinearTensorInterpolationPic2D<T> linearTensorInterpolationPic,
                         const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D> coeff,
                         DenseMatrixT &currentEvaluations,
                         const char* filename, T &L2Error, T &LinftyError)
{
    std::cerr << "Just for testing..." << std::endl;

    int PlotPts1 = currentEvaluations.numRows(), PlotPts2 = currentEvaluations.numCols();
    T h1 = 1./(PlotPts1-1), h2 = 1./(PlotPts2-1);
    std::stringstream PlotFileName;
    PlotFileName << filename << ".dat";
    std::ofstream plotfile(PlotFileName.str().c_str());
    plotfile.precision(6);

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
                T tmp = 0.0;
                T exact= linearTensorInterpolationPic.evaluateInterpolation(x1,x2);

                tmp += coeff * basis.first.generator(xtype1)(x1,j1,k1,0) * basis.second.generator(xtype2)(x2,j2,k2,0);
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

            plotfile << i1 << " " << i2 << " " << appr  << std::endl;
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
plotImageScatterCoeff2D(const LinearTensorInterpolationPic2D<T> linearTensorInterpolationPic,
                        const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D> coeff,
                        const char* filename)
{
    int N1 = linearTensorInterpolationPic.N1;
    int N2 = linearTensorInterpolationPic.N2;

    std::stringstream PlotFileName;
    PlotFileName << filename << ".dat";
    std::ofstream plotfile(PlotFileName.str().c_str());
    plotfile.precision(16);

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


/*
void
plotImageApproximation2D(const LinearTensorInterpolationPic2D<T> linearTensorInterpolationPic,
                         const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D> coeff,
                         T a1, T b1, T a2, T b2, T h, const char* filename)
{
    std::stringstream PlotFileName;
    PlotFileName << filename << ".dat";
    std::ofstream plotfile(PlotFileName.str().c_str());
    plotfile.precision(16);

    for (T x1=a1; x1<=b1; x1+=h) {
        for (T x2=a2; x2<=b2; x2+=h) {
            T appr = 0.0;
            T exact= linearTensorInterpolationPic.evaluateInterpolation(x1,x2);
            for (const_coeff2d_it it = coeff.begin(); it != coeff.end(); ++it) {
                XType xtype1 = (*it).first.index1.xtype;
                XType xtype2 = (*it).first.index2.xtype;
                int j1 = (*it).first.index1.j, j2 = (*it).first.index2.j;
                long k1 = (*it).first.index1.k, k2 = (*it).first.index2.k;

                T coeff = (*it).second;

                appr += coeff * basis.first.generator(xtype1)(x1,j1,k1,0) * basis.second.generator(xtype2)(x2,j2,k2,0);

            }
            appr += linearTensorInterpolationPic.evaluateLiftingFunction(x1,x2);
            exact = std::max(exact,0.); appr = std::max(appr,0.);
            exact = std::min(exact,1.); appr = std::min(appr,1.);
            plotfile << x1 << " " << x2 << " " << exact << " " << appr  << std::endl;
        }
        plotfile << std::endl;
    }
    plotfile.close();
}
*/



/*
 *
    cout << "(0,0): " << linearTensorInterpolPic2D.evaluateLiftingFunction(0.,0.) << " " << linearTensorInterpolPic2D.evaluateInterpolation(0.,0.) << endl;
    cout << "(1,0): " << linearTensorInterpolPic2D.evaluateLiftingFunction(1.,0.) << " " << linearTensorInterpolPic2D.evaluateInterpolation(1.,0.) << endl;
    cout << "(0,1): " << linearTensorInterpolPic2D.evaluateLiftingFunction(0.,1.) << " " << linearTensorInterpolPic2D.evaluateInterpolation(0.,1.) << endl;
    cout << "(1,1): " << linearTensorInterpolPic2D.evaluateLiftingFunction(1.,1.) << " " << linearTensorInterpolPic2D.evaluateInterpolation(1.,1.) << endl;

    cout << "(0.5,0): " << linearTensorInterpolPic2D.evaluateLiftingFunction(0.5,0.)  << " " << linearTensorInterpolPic2D.evaluateInterpolation(0.5,0.) << endl;
    cout << "(1.,0.5): " << linearTensorInterpolPic2D.evaluateLiftingFunction(1.,0.5) << " " << linearTensorInterpolPic2D.evaluateInterpolation(1.,0.5) << endl;
    cout << "(0.5,1.): " << linearTensorInterpolPic2D.evaluateLiftingFunction(0.5,1.) << " " << linearTensorInterpolPic2D.evaluateInterpolation(0.5,1.) << endl;
    cout << "(0,0.5): " << linearTensorInterpolPic2D.evaluateLiftingFunction(0.,0.5)  << " " << linearTensorInterpolPic2D.evaluateInterpolation(0.,0.5) << endl;

    Integral_LinearTensorInterpolPic2D_Test integral_linearTensorInterpolPic2D_test(func_linearTensorInterpolPic2D, basis, basis);
    integral_linearTensorInterpolPic2D_test.quadrature.setOrder(4);
    Coefficients<Lexicographical,T,Index2D> u_test;
    int J = 0;
    for (int k1=basis.mra.rangeI(J).firstIndex(); k1<=basis.mra.rangeI(J).lastIndex(); ++k1) {
        Index1D index1(J,k1,XBSpline);
        for (int k2=basis.mra.rangeI(J).firstIndex(); k2<=basis.mra.rangeI(J).lastIndex(); ++k2) {
            Index1D index2(J,k2,XBSpline);
            Index2D index(index1,index2);
            cout << "... treating index " << index << endl;
            T coeff = integral_linearTensorInterpolPic2D_test(J, k1, XBSpline, 0, J, k2, XBSpline, 0);
            u_test[index] = coeff;
        }
    }



    ofstream plotfile("test.txt");
    for (T x1=0.; x1<=1.; x1+=pow2i<T>(-8)) {
        for (T x2=0.; x2<=1.; x2+=pow2i<T>(-8)) {
            T val = 0.;
            for (const_coeff2d_it it=u.begin(); it!=u.end(); ++it) {
                int k1 = (*it).first.index1.k, k2 = (*it).first.index2.k;
                val += (*it).second * basis.mra.phi(x1,J,k1,0) * basis.mra.phi(x2,J,k2,0);
            }
            plotfile << x1 << " " << x2 << " " << val << " "
                     << linearTensorInterpolPic2D.evaluateInterpolation(x1,x2) << " "
                     << linearTensorInterpolPic2D.evaluateLiftingFunction(x1,x2) << endl;
        }
        plotfile << endl;
    }



    for (T i1=0; i1<PlotPts1; ++i1) {
        T x1 = i1*h1;
        for (T i2=0; i2<PlotPts2; ++i2) {
            T x2 = i2*h2;
            T tmp = 0.0;
            T exact= linearTensorInterpolationPic.evaluateInterpolation(x1,x2);
            for (const_coeff2d_it it = coeff.begin(); it != coeff.end(); ++it) {
                XType xtype1 = (*it).first.index1.xtype;
                XType xtype2 = (*it).first.index2.xtype;
                int j1 = (*it).first.index1.j, j2 = (*it).first.index2.j;
                long k1 = (*it).first.index1.k, k2 = (*it).first.index2.k;

                T coeff = (*it).second;

                tmp += coeff * basis.first.generator(xtype1)(x1,j1,k1,0) * basis.second.generator(xtype2)(x2,j2,k2,0);

            }
            currentEvaluations(i1,i2) += tmp;
            T appr = currentEvaluations(i1,i2) + linearTensorInterpolationPic.evaluateLiftingFunction(x1,x2);

            exact = std::max(exact,0.); appr = std::max(appr,0.);
            exact = std::min(exact,1.); appr = std::min(appr,1.);
            plotfile << i1 << " " << i2 << " " << exact << " " << appr  << std::endl;
        }
        plotfile << std::endl;
    }
    plotfile.close();




    for (T i1=0; i1<PlotPts1; ++i1) {
        T x1 = i1*h1;
        for (T i2=0; i2<PlotPts2; ++i2) {
            T x2 = i2*h2;
            T tmp = 0.0;
            T exact= linearTensorInterpolationPic.evaluateInterpolation(x1,x2);
            for (const_coeff2d_it it = coeff.begin(); it != coeff.end(); ++it) {
                XType xtype1 = (*it).first.index1.xtype;
                XType xtype2 = (*it).first.index2.xtype;
                int j1 = (*it).first.index1.j, j2 = (*it).first.index2.j;
                long k1 = (*it).first.index1.k, k2 = (*it).first.index2.k;

                T coeff = (*it).second;

                tmp += coeff * basis.first.generator(xtype1)(x1,j1,k1,0) * basis.second.generator(xtype2)(x2,j2,k2,0);

            }
            currentEvaluations(i1,i2) += tmp;
            T appr = currentEvaluations(i1,i2) + linearTensorInterpolationPic.evaluateLiftingFunction(x1,x2);

            exact = std::max(exact,0.); appr = std::max(appr,0.);
            exact = std::min(exact,1.); appr = std::min(appr,1.);
            plotfile << i1 << " " << i2 << " " << exact << " " << appr  << std::endl;
        }
        plotfile << std::endl;
    }
    plotfile.close();
 *
 */

