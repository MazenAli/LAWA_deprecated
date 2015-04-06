#include <stdio.h>
#include <cassert>

#include <flens/flens.cxx>

#include <htucker/htucker.h>

#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/constructions/basis.h>
#include <lawa/constructions/interval/multi/basis.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>

typedef double                                                  T;
typedef flens::GeMatrix<flens::FullStorage<T,
                        flens::ColMajor> >                      GeMatrix;
typedef lawa::Basis<double, lawa::Orthogonal, lawa::Interval,
            lawa::Multi>                                        MBasis;
typedef lawa::Coefficients<lawa::Lexicographical,
                            double, lawa::Index1D>              Coefficients;

MBasis   basis(4);
htucker::DimensionIndex di(3, 1);

class HTCoeffFrame
{
    private:
        Coefficients*   U;
        const int       numCols_;
    public:
        typedef typename Coefficients::const_iterator   const_iterator;
        typedef typename Coefficients::iterator         iterator;

        HTCoeffFrame(const int& _numCols_)
                    :numCols_(_numCols_)
        {
            assert(numCols_>0);
            U = new Coefficients [numCols_];
        }

        HTCoeffFrame(const HTCoeffFrame& copy):numCols_(copy.numCols())
        {
            assert(numCols_>0);
            U = new Coefficients [numCols_];
            for (int i=0; i<numCols_; ++i) {
                U[i] = copy[i+1];
            }
        }

        ~HTCoeffFrame()
        {
            if (U)
            {
                delete[] U;
            }
        }

        int
        numCols() const
        {
            return numCols_;
        }

        Coefficients&
        operator[] (const int& col_num)
        {
            assert(col_num>=1 && col_num<=numCols_);
            return U[col_num-1];
        }

        const Coefficients&
        operator[] (const int& col_num) const
        {
            assert(col_num>=1 && col_num<=numCols_);
            return U[col_num-1];
        }

        T&
        operator() (const lawa::Index1D& lambda, const int& col_num)
        {
            assert(col_num>=1 && col_num<=numCols_);
            return U[col_num-1][lambda];
        }
};


class HTCoeffNode : private htucker::HTuckerTreeNode<T>
{
    private:
        const MBasis&                   basis;
        lawa::IndexSet<lawa::Index1D>*  activex;
        int                             numCols;

    public:
        HTCoeffNode(const MBasis& _basis,
                    const htucker::DimensionIndex& _index)
                    :basis(_basis),
                    activex(0),
                    htucker::HTuckerTreeNode<T>(_index),
                    numCols(getUorB().numCols()){}

        HTCoeffNode(const HTCoeffNode& copy):   basis(copy.getBasis()),
                                                numCols(copy.getRank())
        {
            assert(numCols>0);
            activex = new lawa::IndexSet<lawa::Index1D> [numCols];
            for (int i=0; i<numCols; ++i) {
                activex[i] = copy.getActivex(i+1);
            }
        }


        ~HTCoeffNode()
        {
            if (activex) {
                delete[] activex;
            }
        }

        const lawa::IndexSet<lawa::Index1D>&
        getActivex(const int& col_num) const
        {
            assert(col_num>=1 && col_num<=(numCols));
            return activex[col_num-1];
        }

        int
        getRank() const
        {
            return numCols;
        }

        const MBasis&
        getBasis() const
        {
            return basis;
        }

        int
        mapCoeff(const lawa::Index1D& lambda) const
        {
            int j               = lambda.j;
            int k               = lambda.k;
            lawa::XType type    = lambda.xtype;
            int j0              = basis.j0;
            int offset          = 0;

            if (type==lawa::XBSpline) {
                assert(j==j0);
                assert( k>=basis.mra.rangeI(j).firstIndex() &&
                        k<=basis.mra.rangeI(j).lastIndex());
                offset -= basis.mra.rangeI(j).firstIndex();
                ++offset;
                return k+offset;
            } else {
                assert(j>=j0);
                assert(type==lawa::XWavelet);
                assert( k>=basis.rangeJ(j).firstIndex() &&
                        k<=basis.rangeJ(j).lastIndex());
                offset += basis.mra.cardI(j);
                offset -= basis.rangeJ(j).firstIndex();
                ++offset;
                return k+offset;
            }
        }

        void
        addCoeff(const HTCoeffFrame& u)
        {
            using flens::_;
            typedef typename HTCoeffFrame::const_iterator   const_it;

            int numRows = getUorB().numRows();
            int current(0), max(0);

            assert(numCols>0);
            for (int i=1; i<=numCols; ++i) {
                for (const_it it=u[i].cbegin(); it!=u[i].cend(); ++it) {
                    activex[i-1].insert((*it).first);
                    current = mapCoeff((*it).first);
                    if (current>max) {
                        max = current;
                    }
                }
            }

            assert(numRows>0);
            GeMatrix U;
            if (max>numRows) {
                U.resize(max, numCols);
                U(_(1,numRows), _) = getUorB();
            } else {
                U = getUorB();
            }

            for (int i=1; i<=numCols; ++i) {
                for (const_it it=u[i].cbegin(); it!=u[i].cend(); ++it) {
                    current = mapCoeff((*it).first);
                    U(current, i) += (*it).second;
                }
            }

            setUorB(U);
        }

        void
        setCoeff(const HTCoeffFrame& u)
        {
            typedef typename HTCoeffFrame::const_iterator   const_it;

            numCols = u.numCols();
            assert (numCols>0);
            if (activex) {
                delete[] activex;
            }
            activex = new lawa::IndexSet<lawa::Index1D> [numCols];

            int numRows(0), current(0);
            for (int i=1; i<=numCols; ++i) {
                for (const_it it=u[i].cbegin(); it!=u[i].cend(); ++it) {
                    activex[i-1].insert((*it).first);
                    current = mapCoeff((*it).first);
                    if (current>numRows) {
                        numRows = current;
                    }
                }
            }

            assert(numRows>0);
            GeMatrix U(numRows, numCols);
            for (int i=1; i<=numCols; ++i) {
                for (const_it it=u[i].cbegin(); it!=u[i].cend(); ++it) {
                    current = mapCoeff((*it).first);
                    U(current, i) = (*it).second;
                }
            }
            setUorB(U);
        }

        const GeMatrix&
        getFrame() const
        {
            return getUorB();
        }

        void
        printActive()
        {
            std::cout   << "-----   Printing active wavelet indices -----"
                        << std::endl;
            for (int i=0; i<numCols; ++i)
            {
                std::cout << "Column number " << i+1 << std::endl;
                std::cout << activex[i] << std::endl;
            }
            std::cout   << "-----               Done                -----"
                        << std::endl;
        }

};


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
    HTCoeffNode ht(basis, di);
    printf("\n\n---------Testing Map upto j=%d----------------\n\n", j);

    int j0 = basis.j0;
    for (int k=basis.mra.rangeI(j0).firstIndex();
             k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
        lawa::XType type = lawa::XBSpline;
        lawa::Index1D lambda(j0, k, type);
        int coeff = ht.mapCoeff(lambda);
        printf("Lambda is (j,k,xtype) = (%d,%d,%d)\n", j0, k, (int)type);
        printf("Mapped number is coeff = %d\n", coeff);
    }

    printf("\n-------------\n");
    for (int i=basis.j0; i<=j; ++i) {
        for (int k=basis.rangeJ(i).firstIndex();
                 k<=basis.rangeJ(i).lastIndex(); ++k) {
            lawa::XType type = lawa::XWavelet;
            lawa::Index1D lambda(i, k, type);
            int coeff = ht.mapCoeff(lambda);
            printf("Lambda is (j,k,xtype) = (%d,%d,%d)\n", i, k, (int)type);
            printf("Mapped number is coeff = %d\n", coeff);
        }
    }
}


int
main()
{
    const int rank = 3;
    HTCoeffNode U(basis, di);
    HTCoeffFrame u (rank);
    lawa::Index1D       l11(0, 2, lawa::XBSpline),
                        l21(1, 2, lawa::XWavelet),
                        l31(1, 3, lawa::XWavelet),
                        l41(2, 1, lawa::XWavelet);
    lawa::Index1D       l12(0, 1, lawa::XBSpline),
                        l22(1, 1, lawa::XWavelet);
    lawa::Index1D       l13(0, 3, lawa::XBSpline),
                        l23(1, 1, lawa::XWavelet),
                        l33(1, 4, lawa::XWavelet);
    lawa::Index1D       l32(0, 2, lawa::XWavelet),
                        l43(3, 1, lawa::XWavelet);
    u(l11,1) = -1;
    u(l21,1) = 2;
    u(l31,1) = 0.5;
    u(l41,1) = 0.3;
    u(l12,2) = -1.5;
    u(l22,2) = 3;
    u(l13,3) = 1.7;
    u(l23,3) = -0.4;
    u(l33,3) = 7;

    U.setCoeff(u);

    std::cout << U.getFrame() << std::endl;

    u(l32,2) = 0.333;
    u(l43,3) = -0.353;

    U.addCoeff(u);

    std::cout << "After addCoeff" << std::endl;
    std::cout << U.getFrame() << std::endl;

    U.printActive();

//    printBasis(basis, 3);
//    testMap(3);
    return 0;
}
