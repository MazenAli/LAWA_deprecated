#ifndef LAWA_BOX_TENSORBASIS_H
#define LAWA_BOX_TENSORBASIS_H 1

namespace lawa{
    
template<typename FirstBasis, typename SecondBasis>
struct TensorBasis
{
    typedef FirstBasis FirstBasisType;
    typedef SecondBasis SecondBasisType;
    
    TensorBasis(const FirstBasis &_basis1, const SecondBasis &_basis2);

    const FirstBasis &first;
    const SecondBasis &second;
};
    
    
} // namespace lawa

#include <lawa/box/tensorbasis.tcc>


#endif // LAWA_BOX_TENSORBASIS_H