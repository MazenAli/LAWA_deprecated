#include <stdio.h>
#include <cassert>
#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/constructions/basis.h>
#include <lawa/constructions/interval/multi/basis.h>

typedef lawa::Basis<double, lawa::Orthogonal, lawa::Interval,
            lawa::Multi>                                        MBasis;

MBasis   basis(4);

int
mapCoeff(lawa::Index1D lambda)
{
    int j               = lambda.j;
    int k               = lambda.k;
    lawa::XType type    = lambda.xtype;
    int j0              = basis.j0;
    int offset          = 0;

    if (type==lawa::XBSpline) {
        assert(j==j0);
        offset -= basis.mra.rangeI(j).firstIndex();
        ++offset;
        return k+offset;
    } else {
        assert(type==XWavelet);
        offset += basis.mra.cardI(j);
        offset -= basis.rangeJ(j).firstIndex();
        ++offset;
        return k+offset;
    }
}


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
    printf("\n\n---------Testing Map upto j=%d----------------\n\n", j);

    int j0 = basis.j0;
    for (int k=basis.mra.rangeI(j0).firstIndex();
             k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
        lawa::XType type = lawa::XBSpline;
        lawa::Index1D lambda(j0, k, type);
        int coeff = mapCoeff(lambda);
        printf("Lambda is (j,k,xtype) = (%d,%d,%d)\n", j0, k, (int)type);
        printf("Mapped number is coeff = %d\n", coeff);
    }

    printf("\n-------------\n");
    for (int i=basis.j0; i<=j; ++i) {
        for (int k=basis.rangeJ(i).firstIndex();
                 k<=basis.rangeJ(i).lastIndex(); ++k) {
            lawa::XType type = lawa::XWavelet;
            lawa::Index1D lambda(i, k, type);
            int coeff = mapCoeff(lambda);
            printf("Lambda is (j,k,xtype) = (%d,%d,%d)\n", i, k, (int)type);
            printf("Mapped number is coeff = %d\n", coeff);
        }
    }
}


int
main()
{
    int j = 3;

    printBasis(basis, j+1);
    testMap(j);

    return 0;
}
