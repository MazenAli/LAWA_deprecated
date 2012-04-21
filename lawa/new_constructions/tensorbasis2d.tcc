#ifndef LAWA_CONSTRUCTIONS_TENSORBASIS2D_TCC
#define LAWA_CONSTRUCTIONS_TENSORBASIS2D_TCC 1

namespace lawa{
    
template<MethodType Method, typename FirstBasis, typename SecondBasis>
TensorBasis2D<Method, FirstBasis, SecondBasis>::TensorBasis2D(const FirstBasis &_basis1,
                                                              const SecondBasis &_basis2)
    : first(_basis1), second(_basis2)
{
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_TENSORBASIS2D_TCC
