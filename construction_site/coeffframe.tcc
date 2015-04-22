#include <cassert>
#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>


#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFFRAME_TCC
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFFRAME_TCC 1

namespace lawa
{


template <SortingCriterion S, typename T, typename _Index>
CoeffFrame<S, T, _Index>::CoeffFrame(const int& _numCols_)
                                    :numCols_(_numCols_)
{
    assert(numCols_>0);
    U = new lawa::Coefficients<S, T, _Index> [numCols_];
}


template <SortingCriterion S, typename T, typename _Index>
CoeffFrame<S, T, _Index>::CoeffFrame(const CoeffFrame<S, T, _Index>& copy)
                                    :numCols_(copy.numCols())
{
    assert(numCols_>0);
    U = new Coefficients<S, T, _Index> [numCols_];
    for (int i=0; i<numCols_; ++i) {
        U[i] = copy[i+1];
    }
}


template <SortingCriterion S, typename T, typename _Index>
CoeffFrame<S, T, _Index>::~CoeffFrame()
{
    if (U)
    {
        delete[] U;
    }
}


template <SortingCriterion S, typename T, typename _Index>
int
CoeffFrame<S, T, _Index>::numCols() const
{
    return numCols_;
}


template <SortingCriterion S, typename T, typename _Index>
Coefficients<S, T, _Index>&
CoeffFrame<S, T, _Index>::operator[] (const int& col_num)
{
    assert(col_num>=1 && col_num<=numCols_);
    return U[col_num-1];
}


template <SortingCriterion S, typename T, typename _Index>
const Coefficients<S, T, _Index>&
CoeffFrame<S, T, _Index>::operator[] (const int& col_num) const
{
    assert(col_num>=1 && col_num<=numCols_);
    return U[col_num-1];
}


template <SortingCriterion S, typename T, typename _Index>
T&
CoeffFrame<S, T, _Index>::operator() (const _Index& lambda, const int& col_num)
{
    assert(col_num>=1 && col_num<=numCols_);
    return U[col_num-1][lambda];
}


template <SortingCriterion S, typename T, typename _Index>
T
CoeffFrame<S, T, _Index>::operator() (const _Index& lambda, const int& col_num) const
{
    assert(col_num>=1 && col_num<=numCols_);
    return U[col_num-1][lambda];
}


} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFFRAME_TCC
