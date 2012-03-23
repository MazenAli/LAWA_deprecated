#ifndef LAWA_CONSTRUCTIONS_TENSORBASIS3D_TCC
#define LAWA_CONSTRUCTIONS_TENSORBASIS3D_TCC 1

namespace lawa {

template<MethodType Method, typename FirstBasis, typename SecondBasis, typename ThirdBasis>
TensorBasis3D<Method, FirstBasis, SecondBasis, ThirdBasis>::TensorBasis3D(const FirstBasis &_basis1, const SecondBasis &_basis2,
                                                                  const ThirdBasis &_basis3)
    : first(_basis1), second(_basis2), third(_basis3)
{
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_TENSORBASIS3D_TCC
