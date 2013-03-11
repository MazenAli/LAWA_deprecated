namespace lawa {

template <typename T, typename Index, typename RHSType, size_t PDim>
AffineRhs<T, Index,RHSType,PDim>::
AffineRhs(const ThetaStructure<T,PDim>& _thetas, std::vector<RHSType*>& _rhsvec)
 : FlexibleCompoundRhs<T,Index,RHSType>(_rhsvec), thetas(_thetas){}

} // namespace lawa
