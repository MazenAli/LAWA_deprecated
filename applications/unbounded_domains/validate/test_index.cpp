#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef Basis<T,Primal,R,SparseMulti>                                   SparseMW_Basis1D;
typedef Integral<Gauss,SparseMW_Basis1D,SparseMW_Basis1D>               Int;

struct TestIndex1D
{
    XType xtype;
    short j;
    long int k;

    TestIndex1D(size_t tmp=0);
    TestIndex1D(int j, long int k, XType _xtype);
    TestIndex1D(const TestIndex1D &index);
};

TestIndex1D::TestIndex1D(size_t tmp)
: j(0), k(0), xtype(XBSpline)
{
}

TestIndex1D::TestIndex1D(int _j, long int _k, XType _xtype)
: j(_j), k(_k), xtype(_xtype)
{
}

TestIndex1D::TestIndex1D(const TestIndex1D &index)
: j(index.j), k(index.k), xtype(index.xtype)
{
}

struct TestIndex2D
{
    TestIndex1D index1, index2;
};

typedef Coefficients<Lexicographical,double,Index1D>::const_iterator  coeff_iterator;
typedef std::list<Index1D>::const_iterator                       list_iterator;

int main () {
    cout.precision(16);

    SparseMW_Basis1D basis(4,0);

    Int integral(basis,basis);

    int j = 53 /*59*/;
    long int k = 10808639105689173L /*11529215046068477L*/;
    XType xtype = XWavelet;
    Index1D index1d(j,k,xtype);

    cout << "Support:  " << basis.psi.support(j,k) << endl;
    cout << "Tic size: " << basis.psi.tic(j) << endl;
    cout << "Integral: " << integral(j,k,xtype,0, j,k,xtype,0) << endl;

    IndexSet<Index1D> Lambda, LambdaTilde;
    Lambda.insert(index1d);
    cout << Lambda << endl;
    Lambda = C(Lambda, 0.25, basis);
    cout << "Security zone: " << Lambda << endl;

    LambdaTilde = lambdaTilde1d_PDE(index1d, basis, 1, basis.j0, basis.j0+3);
    cout << "LambdaTilde :  " << Lambda << endl;
/*
    TestIndex1D testindex(3,2,XBSpline);
    cout << "Size of test index: " << sizeof(testindex)  << " bytes." << endl;

    Index1D index(2,3,XBSpline);
    cout << "Size of index: " << sizeof(index)  << " bytes." << endl;

    TestIndex2D testindex2d;
    cout << "Size of 2d-testindex: " << sizeof(testindex2d)  << " bytes." << endl;


    Index2D index2d(index,index);
    cout << "Size of 2d-index: " << sizeof(index2d)  << " bytes." << endl << endl << endl;
*/



    long int l = 2L;
    cout << "Size of long int:    " << sizeof(l) << " bytes." << endl;
    float x = 0.1;
    cout << "Size of float:       " << sizeof(x) << " bytes." << endl;
    double y = 0.1;
    cout << "Size of double:      " << sizeof(y) << " bytes." << endl;
    long double z = 0.1L;
    cout << "Size of long double: " << sizeof(z) << " bytes." << endl;

    int int_k   = 2147483647L;
    int int_kP1 = 2147483646L; //2147483648L;
    cout << "int: k =      " << int_k << ", k+1 = " << int_kP1 << endl;

    long int long_int_k   = 9223372036854775806;
    long int long_int_kP1 = 9223372036854775807;
    cout << "long int: k = " << long_int_k << ", k+1 = " << long_int_kP1 << endl;


    cout.precision(20);
    float float_long_int_k = (float)long_int_k;
    cout << "float : k =   " << float_long_int_k << endl;
    cout << "float : k =   " << float_long_int_k-2 << endl;

    return 0;
}

/*
    Coefficients<Lexicographical,double,Index1D>  coeff;
    std::list<Index1D> List;
    for (int j=-1; j<=20; ++j) {
        for (int k=-10000; k<=10000; ++k) {
            coeff[Index1D(j,k,XWavelet)] = rand();
            List.push_back(Index1D(j,k,XWavelet));
        }
    }

    Timer time;
    time.start();
    long double test = 0.;
    for (coeff_iterator it=coeff.begin(); it!=coeff.end(); ++it) {
        test += (long double) (*it).second;
    }
    time.stop();
    cout << "Direct iterating over hash map took: " << time.elapsed() << " seconds, test = " << test << endl;

    time.start();
    test = 0.;
    for (list_iterator it=List.begin(); it!=List.end(); ++it) {
          test += (long double) coeff[(*it)];
    }
    time.stop();
    cout << "Indirect iterating over hash map took: " << time.elapsed() << " seconds, test = " << test << endl;

    int val = -1;
    cout << val << " " << endl;
    val = val << 1;
    cout << val << " " << endl;
*/


