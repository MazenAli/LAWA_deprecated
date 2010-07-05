namespace lawa{
    
template<typename FirstBasis, typename SecondBasis>
TensorBasis<FirstBasis, SecondBasis>::TensorBasis(const FirstBasis &_basis1, const SecondBasis &_basis2)  
    : first(_basis1), second(_basis2)
{
}
    
} // namespace lawa