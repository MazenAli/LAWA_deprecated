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

typedef Integral2D<FullGridGL, PrimalBasis, PrimalBasis>            Integral_LinearTensorInterpolPic2D;

typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

int main (int argc, char *argv[]) {

    //PrimalBasis basis(2,2,0);
    PrimalBasis basis(2,0);
    basis.enforceBoundaryCondition<DirichletBC>();

    LinearTensorInterpolationPic2D<T> linearTensorInterpolPic2D;
    cout << "Reading image data and computing linear tensor interpolation of image data..." << endl;
    linearTensorInterpolPic2D.readPicture("claudiocanutogrey2.txt");
    cout << "... finished." << endl;

    //cout << "Plotting interpolation of image data..." << endl;
    //linearTensorInterpolPic2D.plotInterpolation("claudiocanutogrey2_interpol.txt", 0.002, 0.002);
    //cout << "... finished." << endl;

    Function2D<T> func_linearTensorInterpolPic2D(LinearTensorInterpolationPic2D<T>::evaluateInterpolation,
                                                 LinearTensorInterpolationPic2D<T>::sing_pts_x,
                                                 LinearTensorInterpolationPic2D<T>::sing_pts_y);

    Integral_LinearTensorInterpolPic2D integral_linearTensorInterpolPic2D(func_linearTensorInterpolPic2D, basis, basis);
    integral_linearTensorInterpolPic2D.quadrature.setOrder(4);

    int J=6;

    Coefficients<Lexicographical,T,Index2D> u;
    for (int k1=basis.mra.rangeI(J).firstIndex(); k1<=basis.mra.rangeI(J).lastIndex(); ++k1) {
        Index1D index1(J,k1,XBSpline);
        for (int k2=basis.mra.rangeI(J).firstIndex(); k2<=basis.mra.rangeI(J).lastIndex(); ++k2) {
            Index1D index2(J,k2,XBSpline);
            Index2D index(index1,index2);
            T coeff = integral_linearTensorInterpolPic2D(J, k1, XBSpline, 0, J, k2, XBSpline, 0);
            u[index] = coeff;
        }
    }

    cout << "Size of u: " << u.size() << endl;
    ofstream plotfile("test.txt");
    for (T x1=0.2; x1<=0.8; x1+=0.002) {
        for (T x2=0.2; x2<=0.8; x2+=0.002) {
            T val = 0.;
            for (const_coeff2d_it it=u.begin(); it!=u.end(); ++it) {
                int k1 = (*it).first.index1.k, k2 = (*it).first.index2.k;
                val += (*it).second * basis.mra.phi(x1,J,k1,0) * basis.mra.phi(x2,J,k2,0);
            }
            plotfile << x1 << " " << x2 << " " << val << " "
                     << linearTensorInterpolPic2D.evaluateInterpolation(x1,x2) << endl;
        }
        plotfile << endl;
    }
    plotfile.close();

    //plotLinearInterplolation(A, basis);

    return 0;

}

