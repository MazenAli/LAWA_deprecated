#include <iostream>
#include <fstream>
#include <lawa/lawa.h>
#include <applications/finance/initialconditions/initialconditions.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;

const OptionTypenD optiontype = BasketPut;
T strike = 1.;
T maturity = 1.;
T weight1 = 0.5, weight2 = 0.5;
OptionParameters2D<T,BasketPut> optionparameters(strike, maturity, weight1, weight2, false);

const ProcessType2D  processtype  = BlackScholes2D;
//T r = 0.04; T sigma1 = 0.3, sigma2 = 0.2, rho = 0.;
//T u11 = 1., u12 = 0., u21 = 0., u22 = 1.;
T r = 0.04; T sigma1 = 0.3, sigma2 = 0.2, rho = 0.3;
T u11 = 0.95171801008793943164, u12 = 0.30697366218334239729, u21 = -0.30697366218334239729, u22 = 0.95171801008793943164;
ProcessParameters2D<T,BlackScholes2D>   processparameters(r, sigma1, sigma2, rho, u11, u12, u21, u22);

typedef Basis<T,Orthogonal,Interval,Multi>                          PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;
typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>             Basis2D;

typedef PayoffIntegral2D<optiontype,processtype,Basis2D>            PayoffIntegral;

int main (int argc, char *argv[]) {

    cout.precision(20);
    if (argc!=4) {
        cout << "Usage: " << argv[0] << " d j0 J" << endl;
        return 0;
    }

    int d   = atoi(argv[1]);
    int j0  = atoi(argv[2]);
    int j  = atoi(argv[3]);

    PrimalBasis       basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis  &refinementbasis = basis.refinementbasis;
    Basis2D basis2d(basis,basis);



    PayoffIntegral payoffintegral(optionparameters, processparameters, basis2d);

    cout << "Check for U " << endl;
    DenseMatrixT U(2,2), tU(2,2), Q(2,2), QtU(2,2), UQtU(2,2);
    U  = u11, u12, u21, u22;
    tU = u11, u21, u12, u22;
    Q  = sigma1*sigma1, rho*sigma1*sigma2, rho*sigma1*sigma2, sigma2*sigma2;

    QtU  = Q(1,1)*tU(1,1)+Q(1,2)*tU(2,1), Q(1,1)*tU(1,2)+Q(1,2)*tU(2,2),
           Q(2,1)*tU(1,1)+Q(2,2)*tU(2,1), Q(2,1)*tU(1,2)+Q(2,2)*tU(2,2);
    UQtU = U(1,1)*QtU(1,1)+U(1,2)*QtU(2,1), U(1,1)*QtU(1,2)+U(1,2)*QtU(2,2),
           U(2,1)*QtU(1,1)+U(2,2)*QtU(2,1), U(2,1)*QtU(1,2)+U(2,2)*QtU(2,2);
    cout << "U Q U^T " << UQtU << endl;

    /*
    for (T y1 = -0.5; y1<=0.1; y1+=0.1) {
        T y2 = payoffintegral.adapquad.find_intersectionpoint_y2_given_y1(y1);
        cout << endl;
    }
    for (T y2 = -0.5; y2<=0.1; y2+=0.1) {
        T y1 = payoffintegral.adapquad.find_intersectionpoint_y1_given_y2(y2);
        cout << endl;
    }

    */


    Option2D<T,BasketPut>         basketputoption2d(optionparameters);
    TruncatedBasketPutOption2D<T> truncatedbasketputoption2d;
    truncatedbasketputoption2d.setOption(basketputoption2d);
    truncatedbasketputoption2d.setTransformation(u11, u21, u12, u22);
    truncatedbasketputoption2d.setTruncation(-3., 3., -3., 3., 0, 0.5);

    ofstream plotfile("basketput.dat");
    for (T x1=-3.; x1<=3.; x1+=0.0625) {
        for (T x2=-3.; x2<=3.; x2+=0.0625) {
            plotfile << x1 << " " << x2 << " " << payoffintegral.payoff(x1,x2) << " "
                     << truncatedbasketputoption2d.g_trunc(x1,x2) << endl;
        }
        plotfile << endl;
    }

    Index1D index1(1,2,XBSpline);
    Index1D index2(1,1,XBSpline);
    Index2D index(index1,index2);

    cout << payoffintegral(index) << endl;


    return 0;
}

