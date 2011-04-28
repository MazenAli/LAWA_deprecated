#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEOPERATOR2D_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEOPERATOR2D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa {
    
template <typename T, typename Basis>
struct AdaptiveOperator2D {
               
    virtual T
    operator()(const Index2D &row_index, const Index2D &col_index) = 0;
    
    //CompressionPDE2D<T, Basis> compression;
    
};
    
} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEOPERATOR2D_H