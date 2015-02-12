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
typedef flens::SparseGeMatrix<flens::extensions::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

void
test(std::vector<SparseMatrixT> &vec) {
    cout << "test started..." << endl;
    SparseMatrixT A1(2,2);
    cout << "constructor worked" << endl;
    A1(1,1) = 2.; A1(2,2) =3;
    cout << "assignement worked" << endl;
    A1.finalize();
    cout << "finalized worked" << endl;
    vec.push_back(A1);
    cout << "push back worked" << endl;

    SparseMatrixT A2(2,2);
    A2(1,1) = 1.; A2(2,2) =2;
    A2.finalize();
    vec.push_back(A2);
}

T
test2(const SparseMatrixT &B_store) {
    int N=B_store.numRows();
    cout << "N = " << N << endl;
    SparseMatrixT B1(N,N,B_store.initializer());
    DenseVectorT x(N), y(N);
    for (int i=1; i<=N; ++i) {
        x(i) = 1;
    }
    y = B1*x;
    T ret=0.;
    for (int i=1; i<=N; ++i) {
        ret += y(i);
    }
    return ret;
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

    cout << "Standard constructor for sparse matrices..."  << endl;
    SparseMatrixT A_store(3,3);
    A_store(1,1) = 1.; A_store(2,2) = 1.; A_store(3,3) = 1; A_store(1,2) = 3.;
    //cout << "Size of A_store coordinates: " << A_store._initializer->_coordinates.size() << endl;
    cout << "... worked."  << endl;

    cout << "New constructor for sparse matrices..."  << endl;
    SparseMatrixT A1(3,3,A_store.initializer());
    cout << "... worked."  << endl;

    A1.finalize();
    DenseMatrixT A_dense1;
    densify(cxxblas::NoTrans,A1,A_dense1);
    cout << A_dense1 << endl;

    A_store(1,3) = 4.;
    //A2.initWith(A_store);
    SparseMatrixT A2(3,3,A_store.initializer());
    A2.finalize();
    DenseMatrixT A_dense2;
    densify(cxxblas::NoTrans,A2,A_dense2);
    cout << A_dense2 << endl;

    int N = 5000000;


    SparseMatrixT B_store(N,N);

    for (int i=1; i<=N; ++i) {
        B_store(i,i) = 1.;
    }
    cout << "Starting test2. Hit enter..." << endl;
    getchar();
    cout << test2(B_store) << endl;
    cout << "... finished" << endl;
    getchar();
    cout << "Starting test2. Hit enter..." << endl;
    getchar();
    cout << test2(B_store) << endl;
    cout << "... finished" << endl;
    getchar();


    Timer time;
    time.start();
    SparseMatrixT B1(N,N,B_store.initializer());
    time.stop();
    cout << "Required time for init: " << time.elapsed() << endl;
    getchar();

    B_store.finalize();


    /*
    std::vector<T> ret;
    for (int i=1; i<=N; ++i) {
        ret.push_back(1.);
    }
    getchar();
    std::vector<T> &retPtr = ret;
    cout << retPtr[1] << " " << retPtr[N-1] << endl;
    //cout << (*retPtr)[1] << " " << (*retPtr)[N-1] << endl;
    ret[1] = 2.; ret[N-1] = 3.;
    //cout << (*retPtr)[1] << " " << (*retPtr)[N-1] << endl;
    cout << retPtr[1] << " " << retPtr[N-1] << endl;
    getchar();
    */

    return 0;

}
