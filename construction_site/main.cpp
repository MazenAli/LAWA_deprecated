#include <stdio.h>
#include <flens/flens.cxx>
#include <htucker/htucker.h>
#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/constructions/basis.h>
#include <lawa/constructions/interval/multi/basis.h>
#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptivehelmholtzoperatoroptimized1d.h>
#include <construction_site/coeffframe.h>
#include <construction_site/htcoeffnode.h>

typedef double                                                      T;
typedef lawa::Basis<double, lawa::Orthogonal, lawa::Interval,
            lawa::Multi>                                            MBasis;
typedef lawa::HTCoeffNode<T, MBasis, lawa::Index1D>                 HTCoeffNode;
typedef lawa::CoeffFrame<lawa::Lexicographical, T, lawa::Index1D>   HTCoeffFrame;
typedef lawa::AdaptiveHelmholtzOperatorOptimized1D<T,
        lawa::Orthogonal, lawa::Interval, lawa::Multi>              HeHoOp;

MBasis   basis(4);
htucker::DimensionIndex di(3, 1);


void
printBasis(const MBasis& _basis, int j)
{
    printf("Printing Basis upto level j = %d\n\n", j);
    int j0 = _basis.j0;
    int mracard, card;
    int kl, kr, mkl, mkr;

    for (int i=j0; i<=j; ++i) {
        mracard = _basis.mra.cardI(i);
        card    = _basis.cardJ(i);
        kl      = _basis.rangeJ(i).firstIndex();
        kr      = _basis.rangeJ(i).lastIndex();
        mkl     = _basis.mra.rangeI(i).firstIndex();
        mkr     = _basis.mra.rangeI(i).lastIndex();

        printf("On level j = %d, card = %d, mracard = %d\n", i, card, mracard);
        printf("The range is between %d - %d\n", kl, kr);
        printf("The MRA range is between %d - %d\n", mkl, mkr);
    }
}

void
testMap(int j)
{
    htucker::HTuckerTreeNode<T> htnode(di);
    HTCoeffNode ht(htnode, basis);
    printf("\n\n---------Testing Map upto j=%d----------------\n\n", j);

    int j0 = basis.j0;
    for (int k=basis.mra.rangeI(j0).firstIndex();
             k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
        lawa::XType type = lawa::XBSpline;
        lawa::Index1D lambda(j0, k, type);
        int coeff = lawa::mapCoeff(lambda, basis);
        printf("Lambda is (j,k,xtype) = (%d,%d,%d)\n", j0, k, (int)type);
        printf("Mapped number is coeff = %d\n", coeff);
    }

    printf("\n-------------\n");
    for (int i=basis.j0; i<=j; ++i) {
        for (int k=basis.rangeJ(i).firstIndex();
                 k<=basis.rangeJ(i).lastIndex(); ++k) {
            lawa::XType type = lawa::XWavelet;
            lawa::Index1D lambda(i, k, type);
            int coeff = lawa::mapCoeff(lambda, basis);
            printf("Lambda is (j,k,xtype) = (%d,%d,%d)\n", i, k, (int)type);
            printf("Mapped number is coeff = %d\n", coeff);
        }
    }
}


void
testIterators(HTCoeffNode& ht)
{
    printf("\n\n ------------------- Testing Iterators -------------------\n\n");

    typedef typename HTCoeffNode::const_iterator    const_it;
    typedef typename HTCoeffNode::iterator          iter;
    int numCols = ht.getRank();

    for (int j=1; j<= numCols; ++j) {
        printf("Column j=%d\n", j);
        for (iter it=ht.begin(j); it!=ht.end(j); ++it) {
            std::cout   << "u(["<< (*it) << "], " << j
                        << ") = " << ht((*it), j) << std::endl;
        }
    }

    printf("\n -------------- DONE --------------\n");
}


void
testHTCoeffNode(int j, int rank)
{
    htucker::HTuckerTreeNode<T> htnode(di);
    HTCoeffNode ht(htnode, basis);
    HTCoeffFrame u(rank);
    printf("\n\n---------Testing HTCoeffNode upto j=%d----------------\n", j);
    printf("---------           rank=%d             ----------------\n\n", rank);

    int count(0);
    int j0 = basis.j0;
    for (int l=1; l<=rank; ++l) {
        for (int k=basis.mra.rangeI(j0).firstIndex();
                 k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
            lawa::XType type = lawa::XBSpline;
            lawa::Index1D lambda(j0, k, type);
            ++count;
            u(lambda, l) = count;
        }
        printf("Setting coefficients for scaling fncs\n");
        ht.setCoeff(u);
        std::cout << ht.getFrame() << std::endl;

        printf("\n...Printing Active Coefficients...\n");
        ht.printActive();

        printf("\n-------------\n");
        for (int i=basis.j0; i<=j; ++i) {
            for (int k=basis.rangeJ(i).firstIndex();
                     k<=basis.rangeJ(i).lastIndex(); ++k) {
                lawa::XType type = lawa::XWavelet;
                lawa::Index1D lambda(i, k, type);
                ++count;
                u(lambda, l) = count;
            }
        }
        printf("Adding wavelet coefficients\n");
        ht.addCoeff(u);
        std::cout << ht.getFrame() << std::endl;

        printf("\n...Printing Active Coefficients...\n");
        ht.printActive();
    }
    testIterators(ht);
}

void
apply(HTCoeffNode& u, const T& eps)
{
    T c = 0;
    HeHoOp op(basis, c);

    int numCols = u.getRank();
    HTCoeffFrame U(numCols), ret(numCols);

    typedef typename HTCoeffNode::const_iterator const_it;

    for (int j=1; j<=numCols; ++j) {
        for (const_it it=u.cbegin(j); it!=u.cend(j); ++it) {
            U((*it), j) = u((*it), j);
        }
        op.apply(U[j], eps, ret[j]);
    }

    u.setCoeff(ret);
    printf("Result frame is\n");
    for (int l=1; l<=numCols; ++l) {
        printf("Column %d\n", l);
        std::cout << ret[l] << std::endl;
    }
}


void
testApply(int j, int rank)
{
    htucker::HTuckerTreeNode<T> htnode(di);
    HTCoeffNode ht(htnode, basis);
    HTCoeffFrame u(rank+1);
    printf("\n\n---------Testing apply upto j=%d----------------\n", j);
    printf("---------           rank=%d             ----------------\n\n", rank);

    int count(0);
    T   scale(30);
    int j0 = basis.j0;
    for (int l=1; l<=rank; ++l) {
        for (int k=basis.mra.rangeI(j0).firstIndex();
                 k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
            lawa::XType type = lawa::XBSpline;
            lawa::Index1D lambda(j0, k, type);
            ++count;
            u(lambda, l) = count/scale;
        }
        for (int i=basis.j0; i<=j; ++i) {
            for (int k=basis.rangeJ(i).firstIndex();
                     k<=basis.rangeJ(i).lastIndex(); ++k) {
                lawa::XType type = lawa::XWavelet;
                lawa::Index1D lambda(i, k, type);
                ++count;
                u(lambda, l) = count/scale;
            }
        }
    }

    ht.setCoeff(u);

    printf("Before apply()\n");
    std::cout << ht.getFrame() << std::endl;
    T eps = 1;
    apply(ht, eps);
    printf("After apply()\n");
    std::cout << ht.getFrame() << std::endl;
}


int
main()
{
    //printBasis(basis, 0);
    //testHTCoeffNode(0, 3);
    testApply(0, 3);
    return 0;
}
