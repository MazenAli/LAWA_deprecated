namespace lawa {

template<typename T,size_t PDim>
ThetaStructure<T,PDim>::ThetaStructure()
{}

template<typename T,size_t PDim>
ThetaStructure<T,PDim>::ThetaStructure(const std::vector<ThetaFct>& _thetas)
 : thetas(_thetas)
{}


template<typename T,size_t PDim>
unsigned int
ThetaStructure<T,PDim>::size()
{
	return thetas.size();
}


} // namespace lawa
