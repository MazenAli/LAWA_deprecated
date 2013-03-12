namespace lawa {

template <typename T, typename Index, typename RHSType, size_t PDim>
AffineRhs<T, Index,RHSType,PDim>::
AffineRhs(ThetaStructure<T,PDim>& _thetas, std::vector<RHSType*>& _rhsvec)
 : FlexibleCompoundRhs<T,Index,RHSType>(_rhsvec), thetas(_thetas)
{
	assert(thetas.size() == this->rhsvec.size());
}

template <typename T, typename Index, typename RHSType, size_t PDim>
void
AffineRhs<T, Index,RHSType,PDim>::
set_param(std::array<T,PDim>& mu)
{
	thetas.set_param(mu);
}


template <typename T, typename Index, typename RHSType, size_t PDim>
T
AffineRhs<T, Index,RHSType,PDim>::
operator()(const Index &index)
{
    T val = 0.;

    for(size_t i = 0; i < thetas.size(); ++i){
    	val += thetas.eval(i) * (*(this->rhsvec)[i])(index);
    }

    return val;
}

template <typename T, typename Index, typename RHSType, size_t PDim>
Coefficients<Lexicographical,T,Index>
AffineRhs<T, Index,RHSType,PDim>::
operator()(const IndexSet<Index> &indexset)
{
    Coefficients<Lexicographical,T,Index> r;

    T val;
    for(auto& index : indexset){
    	val = 0.;
        for(size_t i = 0; i < thetas.size(); ++i){
        	val += thetas.eval(i) * (*(this->rhsvec)[i])(index);
        }
        r.insert(std::make_pair(index,val));
    }

    return r;
}

} // namespace lawa
