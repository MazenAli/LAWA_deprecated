#include <iostream>
#include <flens/flens.cxx>
#include <htucker/htucker.h>
#include <lawa/lawa.h>

typedef double                                                      T;
typedef flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> >    DenseM;
typedef flens::DenseVector<flens::Array<T> >                        DenseV;
typedef htucker::HTuckerTreeNode<T>                                 HTNode;
const   flens::Underscore<DenseM::IndexType>                        __;

int
main()
{
    DenseM                  A(3, 3);
    DenseV::View            v = A(__, 1);
    htucker::DimensionIndex j(1);
    j.setValue(1, 3);
    HTNode                  U(j);

    A = 1, 2, 3,
        4, 5, 6,
        7, 8, 9;
    U.setUorB(A);
    DenseM::ConstView       B = U.getUorB();

    std::cout << "B = " << B << std::endl;
    std::cout << "A =" << A << std::endl;

    return 0;
}
