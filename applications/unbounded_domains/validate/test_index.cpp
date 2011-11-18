#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

struct TestIndex1D
{
    XType xtype;
    short j;
    long int k;
    size_t * const p_hash_value;

    TestIndex1D(size_t tmp=0);
    TestIndex1D(int j, long int k, XType _xtype, size_t tmp = 0);
    TestIndex1D(const TestIndex1D &index, size_t tmp=0);
};

TestIndex1D::TestIndex1D(size_t tmp)
: j(0), k(0), xtype(XBSpline), p_hash_value(&tmp)
{
    cout << "Adress of p_hash_value: " << p_hash_value << endl;
}

TestIndex1D::TestIndex1D(int _j, long int _k, XType _xtype, size_t tmp)
: j(_j), k(_k), xtype(_xtype), p_hash_value(&tmp)
{
    cout << "Adress of p_hash_value: " << p_hash_value << endl;
}

TestIndex1D::TestIndex1D(const TestIndex1D &index, size_t tmp)
: j(index.j), k(index.k), xtype(index.xtype), p_hash_value(&tmp)
{
    cout << "Adress of p_hash_value: " << p_hash_value << endl;
}


typedef Coefficients<Lexicographical,double,Index1D>::const_iterator  coeff_iterator;
typedef std::list<Index1D>::const_iterator                       list_iterator;

int main () {

    TestIndex1D testindex;
    cout << "hash_value:  " << *(testindex.p_hash_value) << endl;
    *(testindex.p_hash_value) = (size_t)1;
    cout << "hash_value:  " << *(testindex.p_hash_value) << endl << endl;

    TestIndex1D testindex2(3,2,XBSpline);
    cout << "hash_value:  " << *(testindex.p_hash_value) << endl;
    *(testindex.p_hash_value) = (size_t)2;
    cout << "hash_value:  " << *(testindex.p_hash_value) << endl << endl;

    TestIndex1D testindex3(testindex2);
    cout << "hash_value:  " << *(testindex.p_hash_value) << endl;
    *(testindex.p_hash_value) = (size_t)3;
    cout << "hash_value:  " << *(testindex.p_hash_value) << endl << endl;
/*
    Index1D index(2,3,XBSpline);
    cout << "Size of index: " << sizeof(index)  << " bytes." << endl;

    Index2D index2d(index,index);
    cout << "Size of 2d-index: " << sizeof(index2d)  << " bytes." << endl << endl << endl;


    long int l = 2L;
    cout << "Size of long int:    " << sizeof(l) << " bytes." << endl;
    float x = 0.1;
    cout << "Size of float:       " << sizeof(x) << " bytes." << endl;
    double y = 0.1;
    cout << "Size of double:      " << sizeof(y) << " bytes." << endl;
    long double z = 0.1L;
    cout << "Size of long double: " << sizeof(z) << " bytes." << endl;

    int int_k   = 2147483647L;
    int int_kP1 = 2147483648L;
    cout << "int: k =      " << int_k << ", k+1 = " << int_kP1 << endl;

    long int long_int_k   = 9223372036854775806;
    long int long_int_kP1 = 9223372036854775807;
    cout << "long int: k = " << long_int_k << ", k+1 = " << long_int_kP1 << endl;


    cout.precision(20);
    float float_long_int_k = (float)long_int_k;
    cout << "float : k =   " << float_long_int_k << endl;
    cout << "float : k =   " << float_long_int_k-2 << endl;
*/

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


