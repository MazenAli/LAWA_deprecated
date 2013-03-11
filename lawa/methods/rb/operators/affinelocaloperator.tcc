namespace lawa {

template <typename Index, typename LocalOperatorType, size_t PDim>
AffineLocalOperator<Index,LocalOperatorType,PDim>::
AffineLocalOperator(const ThetaStructure<T,PDim>& _thetas, std::vector<LocalOperatorType*>& _localops)
 : FlexibleCompoundLocalOperator<Index,LocalOperatorType>(_localops), thetas(_thetas){}

} // namespace lawa
