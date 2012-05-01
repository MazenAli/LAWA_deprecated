#include <iostream>
#include <lawa/lawa.h>
#include <tr1/unordered_set>

using namespace std;
using namespace lawa;

//#define P 6291469
#define P 111

struct index2d_hashfunction1
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

struct index2d_hashfunction2
{
    inline
    size_t operator()(const Index2D& index) const
    {
        size_t l1 = (1L << (index.index1.j+index.index1.xtype) ) + index.index1.k;
        size_t l2 = (1L << (index.index2.j+index.index2.xtype) ) + index.index2.k;
        size_t s1 = l1;
        size_t s2 = l1+l2;
        size_t TwoP = 2*P;
        return (((((s2+1)%(TwoP)) * (s2 % TwoP)) % TwoP)/2 + s1 % P) % P;
    }
};

struct index3d_hashfunction1
{
    inline
    size_t operator()(const Index3D& index) const
    {
        size_t l1 = (1L << (index.index1.j+index.index1.xtype) ) + index.index1.k;
        size_t l2 = (1L << (index.index2.j+index.index2.xtype) ) + index.index2.k;
        size_t l3 = (1L << (index.index3.j+index.index3.xtype) ) + index.index3.k;
        size_t s1 = l1;
        size_t s2 = l1+l2;
        size_t s3 = l1+l2+l3;
        return ( ((s3+2)*(s3+1)*s3)/6 + ((s2+1)*s2)/2 + s1) % P;
    }
};

struct index3d_hashfunction2
{
    inline
    size_t operator()(const Index3D& index) const
    {
        size_t l1 = (1L << (index.index1.j+index.index1.xtype) ) + index.index1.k;
        size_t l2 = (1L << (index.index2.j+index.index2.xtype) ) + index.index2.k;
        size_t l3 = (1L << (index.index3.j+index.index3.xtype) ) + index.index3.k;
        size_t s1 = l1;
        size_t s2 = l1+l2;
        size_t s3 = l1+l2+l3;
        size_t TwoP = 2*P;
        size_t SixP = 6*P;
        return (   ( ( ((s3+2)%SixP)*((s3+1)%SixP)*(s3%SixP) )%SixP )/6
                 + ( ( ((s2+1)%TwoP)*(s2%TwoP) ) % TwoP )/2 + s1%P ) % P;

    }
};

int main (int argc, char *argv[]) {

    int J=10;

/*
    index2d_hashfunction1  index2d_hasher1;
    index2d_hashfunction2  index2d_hasher2;

    for (int j1=0; j1<=J; ++j1) {
        for (int k1=1; k1<=pow2i<int>(j1); ++k1) {
            XType xtype1 = (j1 == 0) ? XBSpline : XWavelet;
            Index1D index1(j1,k1,xtype1);
            for (int j2=0; j2<=J-j1; ++j2) {
                for (int k2=1; k2<=pow2i<int>(j2); ++k2) {
                    XType xtype2 = (j2 == 0) ? XBSpline : XWavelet;
                    Index1D index2(j2,k2,xtype2);
                    Index2D index(index1,index2);
                    size_t val1 = index2d_hasher1.operator()(index);
                    size_t val2 = index2d_hasher2.operator()(index);
                    if (val1 != val2) {
                        cout << "Mismatch for " <<  index << ": " << val1 << " " << val2 << endl;
                    }
                }
            }
        }
    }


    cout << "Hit enter to continue." << endl;
    getchar();
*/
    index3d_hashfunction1  index3d_hasher1;
    index3d_hashfunction2  index3d_hasher2;


    Coefficients<Lexicographical,double,Index3D> v1(255), v2(700);
    IndexSet<Index3D> Lambda;

    cout << "Number of buckets: " << v1.bucket_count() << endl;
    cout << "Number of buckets: " << v2.bucket_count() << endl;

    for (int j1=0; j1<=J; ++j1) {
        for (int k1=1; k1<=pow2i<int>(j1); ++k1) {
            XType xtype1 = (j1 == 0) ? XBSpline : XWavelet;
            Index1D index1(j1,k1,xtype1);
            for (int j2=0; j2<=J-j1; ++j2) {
                for (int k2=1; k2<=pow2i<int>(j2); ++k2) {
                    XType xtype2 = (j2 == 0) ? XBSpline : XWavelet;
                    Index1D index2(j2,k2,xtype2);
                    for (int j3=0; j3<=J-j1-j2; ++j3) {
                        for (int k3=1; k3<=pow2i<int>(j3); ++k3) {
                            XType xtype3 = (j3 == 0) ? XBSpline : XWavelet;
                            Index1D index3(j3,k3,xtype3);
                            Index3D index(index1,index2,index3);
                            size_t val1 = index3d_hasher1.operator()(index);
                            size_t val2 = index3d_hasher2.operator()(index);
                            //v1[index] = 1.;
                            Lambda.insert(index);
                            if (val1 != val2) {
                                cout << index << ": " << val1 << " " << val2 << endl;
                            }
                        }
                    }
                }
            }
        }
    }
    cout << "Size of v: " << v1.size() << ", storage: " << sizeof(v1) << endl;
    cout << "Size of Lambda: " << Lambda.size() << ", storage: " << sizeof(Lambda) << endl;
    getchar();

    return 0;
}
