#include <cstddef>
#include <flens/flens.cxx>
#include <htucker/htucker.h>
#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <construction_site/coeffframe.h>

#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HTCOEFFNODE_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HTCOEFFNODE_H 1

namespace lawa
{

template <typename T, typename _Basis, typename _Index>
class HTCoeffNode
{
    private:
        htucker::HTuckerTreeNode<T>&    htnode;
        const _Basis&                   basis;
        IndexSet<_Index>*               activex;
        std::size_t                     numCols_;

    public:
        typedef typename IndexSet<_Index>::const_iterator   const_iterator;
        typedef typename IndexSet<_Index>::iterator         iterator;
        typedef _Index                                      IndexType;

        // Constructors
        HTCoeffNode(htucker::HTuckerTreeNode<T>& _htnode,
                    const _Basis& _basis);

        HTCoeffNode(const HTCoeffNode<T, _Basis, _Index>& copy);

        ~HTCoeffNode();

        // Get methods
        const htucker::HTuckerTreeNode<T>&
        getNode() const;

        const IndexSet<_Index>&
        getActivex(const std::size_t col_num) const;

        std::size_t
        numCols() const;

        const _Basis&
        getBasis() const;

        const flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >&
        getFrame() const;

        T
        operator() (const _Index& lambda, const std::size_t col_num) const;

        void
        printActive();

        bool
        isActive(const _Index& lambda, const std::size_t col_num) const;

        // Iterators over column
        iterator
        begin(const std::size_t j);

        const_iterator
        cbegin(const std::size_t j) const;

        iterator
        end(const std::size_t j);

        const_iterator
        cend(const std::size_t j) const;

        // Set methods
        template<SortingCriterion S>
        void
        setCoeff(const CoeffFrame<S, T, _Index>& u);

        template<SortingCriterion S>
        void
        addCoeff(const CoeffFrame<S, T, _Index>& u);
};


template <typename _Index, typename _Basis>
unsigned long
mapCoeff(const _Index&, const _Basis&);

template <typename _Basis>
unsigned long
mapCoeff(const Index1D& lambda, const _Basis& basis);


} // namespace lawa

#include <construction_site/htcoeffnode.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HTCOEFFNODE_H
