#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>

#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFFRAME_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFFRAME_H 1

namespace lawa
{


template <SortingCriterion S, typename T, typename _Index>
class CoeffFrame
{
    private:
        Coefficients<S, T, _Index>*   U;
        const int       numCols_;
    public:
        typedef typename Coefficients<S, T, _Index>::const_iterator  const_iterator;
        typedef typename Coefficients<S, T, _Index>::iterator        iterator;
        typedef _Index                                               IndexType;

        CoeffFrame(const int _numCols_);

        CoeffFrame(const CoeffFrame<S, T, _Index>& copy);

        ~CoeffFrame();

        // Get methods
        int
        numCols() const;

        const Coefficients<S, T, _Index>&
        operator[] (const int col_num) const;

        T
        operator() (const _Index& lambda, const int col_num) const;

        // Set methods
        Coefficients<S, T, _Index>&
        operator[] (const int col_num);

        T&
        operator() (const _Index& lambda, const int col_num);
};


} // namespace lawa

#include <construction_site/coeffframe.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFFRAME_H
