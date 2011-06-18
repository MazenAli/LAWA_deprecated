#include <iostream>
#include <fstream>
#include <lawa/lawa.h>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>


using namespace std;
using namespace lawa;
using namespace boost::numeric::ublas;

typedef double T;

typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;


void
test(std::vector<SparseMatrixT> &vec) {
    SparseMatrixT A1(2,2);
    A1(1,1) = 2.; A1(2,2) =3;
    A1.finalize();
    vec.push_back(A1);

    SparseMatrixT A2(2,2);
    A2(1,1) = 1.; A2(2,2) =2;
    A2.finalize();
    vec.push_back(A2);
}


int main () {

    std::vector<SparseMatrixT> vec;
    test(vec);
    DenseVectorT x(2);
    x = 1, 1;
    DenseVectorT A1x(2), A2x(2);
    A1x = vec[0]*x;
    cout << A1x << endl;
    A2x = vec[1]*x;
    cout << A2x << endl;


    std::vector<int*> vec2;
    for (int i=0; i<4; ++i) {
        int *pair1 = new int[2];
        pair1[0] = i; pair1[1] = i+1;
        vec2.push_back(pair1);
    }
    for (int i=0; i<4; ++i) {
        cout << vec2[i][0] << " " << vec2[i][1] << endl;
    }

    mapped_vector<double> v (10000000);
    v.insert_element (101,2.) = 1;
    v.insert_element (106,0.) = 0;

    for (mapped_vector<double>::const_iterator it=v.begin(); it!=v.end(); ++it) {
        std::cout << *it << std::endl;
    }
    for (int i=1000; i<=1000000; ++i) {
        if (v(i)!=0)
        std::cout << "Value at i=" << i  <<": "<< v(i) << std::endl;
    }

    getchar();
}
