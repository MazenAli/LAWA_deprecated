#include <iostream>
#include <lawa/lawa.h>
#include <tr1/unordered_set>

using namespace std;
using namespace lawa;

typedef double T;
typedef Basis<T,Primal,Interval,Dijkema>                            PrimalBasis;

#define P 111

struct index_hashfunction_test
{
    inline
    size_t operator()(const Index2D& index) const
    {
        size_t l1 = (1L << (index.index1.j+index.index1.xtype) ) + index.index1.k;
        size_t l2 = (1L << (index.index2.j+index.index2.xtype) ) + index.index2.k;
        size_t s1 = l1;
        size_t s2 = l1+l2;
        return (((s2+1)*s2)/2 + s1) % P;
    }
};

struct index_hashfunction_test2
{
    index_hashfunction_test2(void) {
        std::cerr << "  Constructor called." << std::endl;
    }

    inline
    size_t operator()(const Index2D& index) const
    {
        size_t l1 = (1L << (index.index1.j+index.index1.xtype) ) + index.index1.k;
        size_t l2 = (1L << (index.index2.j+index.index2.xtype) ) + index.index2.k;
        size_t s1 = l1;
        size_t s2 = l1+l2;
        size_t TwoP = 2*P;
        //return ( (((s2+1)*s2)%TwoP)/2 + s1%P) % P;
        return (((((s2+1)%(TwoP)) * (s2 % TwoP)) % TwoP)/2 + s1 % P) % P;
    }
};

int main (int argc, char *argv[]) {

    int d =2;
    int d_=2;
    int j0=2;
    int J=6;

    PrimalBasis basis(d,d_,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    cout << "Size of size_t: " << sizeof(size_t) << ", size of uns. long int: " << sizeof(unsigned long int) << endl;

    long twoP = 2*P;
    for (long i=1; i<=50; ++i) {
        long val1 = (i / 2) % P;
        long val2 = ( i % (twoP) ) / 2;
        cout << "i = " << i << ": " << val1 << ", " << val2 << endl;
        //if (val1!=val2) cout << "i = " << i << ": " << val1 << ", " << val2 << endl;
    }

    /*
    unsigned long largeint = 4294967295;
    size_t        largeint2 = 4294967295;
    cout << largeint << " " << largeint << endl;
    largeint*=2;
    largeint2*=2;
    cout << largeint << " " << largeint << endl;
    */
    index_hashfunction<Index1D> hasher1d;
    index_hashfunction_test  refhasher2d;
    index_hashfunction_test2 hasher2d;

    Index1D index1(40,1,XWavelet);
    Index1D index2(40,120121020,XWavelet);
    Index2D index(index1,index2);

    cout << index << " " << hasher2d.operator()(index) << endl;


    /*
    for (long k=basis.mra.rangeI(j0).firstIndex(); k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
        Index1D index(j0,k,XBSpline);
        cout << index << " " << hasher1d.operator()(index) << endl;
    }
    for (int j=j0; j<=J; ++j) {
        for (long k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
            Index1D index(j,k,XWavelet);
            cout << index << " " << hasher1d.operator()(index) << endl;
        }
    }
    */
    /*
    std::tr1::unordered_set<unsigned long int> set_of_hashvalues;
    for (long k1=basis.mra.rangeI(j0).firstIndex(); k1<=basis.mra.rangeI(j0).lastIndex(); ++k1) {
        for (long k2=basis.mra.rangeI(j0).firstIndex(); k2<=basis.mra.rangeI(j0).lastIndex(); ++k2) {
            Index1D row(j0,k1,XBSpline);
            Index1D col(j0,k2,XBSpline);
            Index2D index(row,col);
            size_t refhashvalue = refhasher2d.operator()(index);
            size_t hashvalue = hasher2d.operator()(index);
            cout << index << " : " << refhashvalue << " " << hashvalue << endl;
            if (set_of_hashvalues.count(refhashvalue)==0) {
                set_of_hashvalues.insert(refhashvalue);
            }
            else {
                cout << "-> Collision" << endl;
            }

        }
        for (int i2=0; i2<=J; ++i2) {
            for (long k2=basis.rangeJ(j0+i2).firstIndex(); k2<=basis.rangeJ(j0+i2).lastIndex(); ++k2) {
                Index1D row(j0,k1,XBSpline);
                Index1D col(j0+i2,k2,XWavelet);
                Index2D index(row,col);
                size_t refhashvalue = refhasher2d.operator()(index);
                size_t hashvalue = hasher2d.operator()(index);
                cout << index << " : " << refhashvalue << " " << hashvalue << endl;
                if (set_of_hashvalues.count(refhashvalue)==0) {    set_of_hashvalues.insert(refhashvalue); }
                else                                       {    cout << "-> Collision" << endl; }
                Index2D index2(col,row);
                refhashvalue = refhasher2d.operator()(index2);
                hashvalue = hasher2d.operator()(index2);
                cout << index << " : " << refhashvalue << " " << hashvalue << endl;
                if (set_of_hashvalues.count(refhashvalue)==0) {    set_of_hashvalues.insert(refhashvalue); }
                else                                       {    cout << "-> Collision" << endl; }
            }
        }
    }
    for (int i1=0; i1<=J; ++i1) {
        for (long k1=basis.rangeJ(j0+i1).firstIndex(); k1<=basis.rangeJ(j0+i1).lastIndex(); ++k1) {
            for (int i2=0; i1+i2<=J; ++i2) {
                for (long k2=basis.rangeJ(j0+i2).firstIndex(); k2<=basis.rangeJ(j0+i2).lastIndex(); ++k2) {
                    Index1D row(j0+i1,k1,XWavelet);
                    Index1D col(j0+i2,k2,XWavelet);
                    Index2D index(row,col);
                    size_t refhashvalue = refhasher2d.operator()(index);
                    size_t hashvalue = hasher2d.operator()(index);
                    cout << index << " : " << refhashvalue << " " << hashvalue << endl;
                    if (set_of_hashvalues.count(refhashvalue)==0) {    set_of_hashvalues.insert(refhashvalue); }
                    else                                       {    cout << "-> Collision" << endl; }
                }
            }
        }
    }
    */
    return 0;
}
