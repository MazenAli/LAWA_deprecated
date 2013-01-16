#include <iostream>
#include <fstream>
#include <lawa/lawa.h>
#include <applications/finance/initialconditions/initialconditions.h>

using namespace std;
using namespace lawa;

typedef long double T;

typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;

typedef Basis<T,Orthogonal,Interval,Multi>                          PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;
typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>             Basis2D;

//typedef PayoffIntegral2D<FullGridGL,Basis2D,TruncatedBasketPutOption2D<T> > PayoffIntegral;
typedef PayoffIntegral2D<FullGridGL,Basis2D,TruncatedSumOfPutsOption2D<T> > PayoffIntegral;

typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

const OptionTypenD optiontype = SumOfPuts;
T strike = 1.;
T maturity = 1.;
T weight1 = 0.5, weight2 = 0.5;
//OptionParameters2D<T,BasketPut> optionparameters(strike, maturity, weight1, weight2, false);
OptionParameters2D<T,SumOfPuts> optionparameters(strike, strike, maturity, weight1, weight2, false);

const ProcessType2D  processtype  = BlackScholes2D;
//T r = 0.04; T sigma1 = 0.3, sigma2 = 0.2, rho = 0.;
//T u11 = 1., u12 = 0., u21 = 0., u22 = 1.;
T r = 0.04; T sigma1 = 0.3, sigma2 = 0.2, rho = 0.3;
T u11 = 0.95171801008793943164, u12 = 0.30697366218334239729, u21 = -0.30697366218334239729, u22 = 0.95171801008793943164;
T    critical_line_x1 = 0.4;
bool critical_above_x1 = true;

ProcessParameters2D<T,BlackScholes2D>   processparameters(r, sigma1, sigma2, rho, u11, u12, u21, u22);


T
evaluate(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
         Coefficients<Lexicographical,T,Index2D> &v, T x1, T x2);

int main (int argc, char *argv[]) {

    cout.precision(20);
    if (argc!=5) {
        cout << "Usage: " << argv[0] << " d j0 J order" << endl;
        return 0;
    }

    int d     = atoi(argv[1]);
    int j0    = atoi(argv[2]);
    int j     = atoi(argv[3]);
    int order = atoi(argv[4]);

    PrimalBasis       basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis  &refinementbasis = basis.refinementbasis;
    Basis2D basis2d(basis,basis);

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

    T left_x1 = -3., right_x1 = 3.;
    T left_x2 = -3., right_x2 = 3.;

    Option2D<T,optiontype>         option2d(optionparameters);
    //TruncatedBasketPutOption2D<T> truncatedoption2d;
    TruncatedSumOfPutsOption2D<T> truncatedoption2d;
    truncatedoption2d.setOption(option2d);
    //truncatedoption2d.setTransformation(u11, u21, u12, u22);
    truncatedoption2d.setTruncation(left_x1, right_x1, left_x2, right_x2, 0, 0.5, 20.);
    truncatedoption2d.setCriticalLine_x1(critical_line_x1, critical_above_x1);

    PayoffIntegral payoffIntegral(basis2d, truncatedoption2d,
                                  left_x1, right_x1, left_x2, right_x2, false, 1e-1, order);

    Coefficients<Lexicographical,T,Index2D> v;
    getSparseGridVector(basis2d, v, j, (T)0.);

    Timer time;
    time.start();
    for (coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second = payoffIntegral((*it).first);
    }
    time.stop();
    cout.precision(10);
    cout << "v.size() = " << v.size() << ", time elapsed: " << time.elapsed() << endl;

    ofstream plotfile("basketput.dat");
    T maxerror = 0.;
    for (T x1=left_x1; x1<=right_x1; x1+=0.03125) {
        for (T x2=left_x2; x2<=right_x2; x2+=0.03125) {
            T exact = truncatedoption2d.payoff(x1,x2);
            T approx = 0.;
            approx = evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, v, x1, x2);
            plotfile << x1 << " " << x2 << " " << exact << " " << approx << endl;
            maxerror = std::max(maxerror, fabs(approx-exact));
        }
        plotfile << endl;
    }
    cout << "Maximum error: " << maxerror << endl;
    plotfile.close();

    return 0;
}

T
evaluate(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
         Coefficients<Lexicographical,T,Index2D> &v, T x1, T x2)
{
    T RightmLeft_x1 = right_x1-left_x1, SqrtRightmLeft_x1 = std::sqrt(right_x1-left_x1);
    T RightmLeft_x2 = right_x2-left_x2, SqrtRightmLeft_x2 = std::sqrt(right_x2-left_x2);

    T ret = 0.;
    for (const_coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        int   j1 = (*it).first.index1.j,     j2 = (*it).first.index2.j;
        int   k1 = (*it).first.index1.k,     k2 = (*it).first.index2.k;
        XType e1 = (*it).first.index1.xtype, e2 = (*it).first.index2.xtype;

        T val_x1 = (1./SqrtRightmLeft_x1) * basis2d.first.generator(e1).operator()((x1-left_x1)/(RightmLeft_x1),j1,k1,0);
        T val_x2 = (1./SqrtRightmLeft_x2) * basis2d.second.generator(e2).operator()((x2-left_x2)/(RightmLeft_x2),j2,k2,0);

        ret += (*it).second * val_x1 * val_x2;
    }
    return ret;
}

