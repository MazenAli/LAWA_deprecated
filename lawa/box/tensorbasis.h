#ifndef LAWA_BOX_TENSORBASIS_H
#define LAWA_BOX_TENSORBASIS_H 1

namespace lawa{
    
template<typename FirstBasis, typename SecondBasis>
struct TensorBasis
{
    TensorBasis(const FirstBasis &_basis1, SecondBasis &_basis2);

    FirstBasis basis1;
    SecondBasis basis2;
};
    
    
} // namespace lawa

#include <lawa/box/tensorbasis.tcc>


#endif // LAWA_BOX_TENSORBASIS_H