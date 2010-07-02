namespace lawa{
    
template<typename FirstBasis, typename SecondBasis>
TensorBasis<FirstBasis, SecondBasis>::TensorBasis(const FirstBasis &_basis1, const SecondBasis &_basis2, int _J1, int _J2)  
    : basis1(_basis1), basis2(_basis2)
{
}
    
} // namespace lawa