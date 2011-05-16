#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

// Basis definitions
typedef Basis<T,Primal,R,CDF>                   Basis1D;

// Operator definitions
typedef WeightedHelmholtzOperator1D<T, Basis1D>             WeightedPDEOp;
typedef WeightedSobolevMidPointPreconditioner1D<T,Basis1D>  WeightedPreconditioner;
typedef H1NormPreconditioner1D<T,Basis1D>                   Preconditioner;

// FLENS definitions
typedef flens::DenseVector<flens::Array<T> >    DenseVectorT;

const T eta=2.;

int main(int argc, char *argv[]) {

    cout.precision(16);
    if (argc != 5) {
        cout << "usage " << argv[0] << " d d_ j0 k" << endl; exit(1);
    }
    int d =atoi(argv[1]);
    int d_=atoi(argv[2]);
    int j0=atoi(argv[3]);
    int k =atoi(argv[4]);

    Basis1D basis(d,d,j0);

    ExponentialWeightFunction1D<T> expweight;
    Function<T>                    expweightFct(expweight.weight,expweight.singularPoints);
    expweight.setEta(eta);

    WeightedPDEOp           weighted_op(basis,1.,expweightFct,4);
    WeightedPreconditioner  weighted_P(basis,expweightFct,1);
    Preconditioner          P(basis);
    IntegralExpWeight<Gauss,Basis1D,Basis1D> integral_expweight(basis,basis,eta);

    Index1D row1(j0,k,XBSpline), col1(j0,k+1,XBSpline);
    Index1D row2(j0,k,XWavelet), col2(j0,k,XWavelet);

    cout << "Test1: " << endl;
    cout << weighted_P(row1) * weighted_op(row1,col1) * weighted_P(col1) << ": "
         << weighted_P(row1) << " " <<  weighted_op(row1,col1) << " " << weighted_P(col1) << endl;

    T tmp1=0.;
    tmp1 += integral_expweight(row1.j,row1.k,row1.xtype,0,col1.j,col1.k,col1.xtype,0);
    tmp1 += integral_expweight(row1.j,row1.k,row1.xtype,1,col1.j,col1.k,col1.xtype,1);
    cout << P(row1) * tmp1 * P(col1) << ": "
         << P(row1) << " " <<  tmp1 << " " << P(col1)  << endl;

    expweight.setPrecPoints(0.,0.);
    cout << "Test2: " << endl;
    cout << weighted_P(row2) * weighted_op(row2,col2) * weighted_P(col2) << ": "
         << weighted_P(row2) << " " <<  weighted_op(row2,col2) << " " << weighted_P(col2) << endl;

    T tmp2=0.;
    tmp2 += integral_expweight(row2.j,row2.k,row2.xtype,0,col2.j,col2.k,col2.xtype,0);
    tmp2 += integral_expweight(row2.j,row2.k,row2.xtype,1,col2.j,col2.k,col2.xtype,1);
    cout << P(row2) * tmp2 * P(col2) << ": "
         << P(row2) << " " <<  tmp2 << " " << P(col2)  << endl;

    return 0;
}
