namespace lawa {

template <typename Index, typename LocalOperatorType, size_t PDim>
AffineLocalOperator<Index,LocalOperatorType,PDim>::
AffineLocalOperator(const ThetaStructure<T,PDim>& _thetas, std::vector<LocalOperatorType*>& _localops)
 : FlexibleCompoundLocalOperator<Index,LocalOperatorType>(_localops), thetas(_thetas)
{
	assert(thetas.size() == this->localops.size());
}

template <typename Index, typename LocalOperatorType, size_t PDim>
void
AffineLocalOperator<Index,LocalOperatorType,PDim>::eval(const Coefficients<Lexicographical,T,Index> &v,
		 	 	 	 	 	 	 	 	 	 	 	 	 	  Coefficients<Lexicographical,T,Index> &Av)
{
	Coefficients<Lexicographical,typename LocalOperatorType::T,Index> tmp(Av);

    for(size_t i = 0; i < thetas.size(); ++i){
    	tmp.setToZero();
		this->localops[i]->eval(v, tmp);
		Av += thetas.eval(i) * tmp;
	}
}

template <typename Index, typename LocalOperatorType, size_t PDim>
template <typename Preconditioner>
void
AffineLocalOperator<Index,LocalOperatorType,PDim>::eval(Coefficients<Lexicographical,T,Index> &v,
														Coefficients<Lexicographical,T,Index> &Av, Preconditioner &P)
{
	Coefficients<Lexicographical,typename LocalOperatorType::T,Index> tmp(Av);

    for(size_t i = 0; i < thetas.size(); ++i){
    	tmp.setToZero();
    	this->localops[i]->eval(v, tmp, P);
		Av += thetas.eval(i) * tmp;
	}
}

template <typename Index, typename LocalOperatorType, size_t PDim>
template <typename RightPrec, typename LeftPrec>
void
AffineLocalOperator<Index,LocalOperatorType,PDim>::eval(Coefficients<Lexicographical,T,Index> &v,
	 Coefficients<Lexicographical,T,Index> &Av, RightPrec& rightP, LeftPrec& leftP)
 {
	Coefficients<Lexicographical,typename LocalOperatorType::T,Index> tmp(Av);

	 for(size_t i = 0; i < thetas.size(); ++i){
		tmp.setToZero();
		this->localops[i]->eval(v, tmp, rightP, leftP);
		Av += thetas.eval(i) * tmp;
	}
 }

} // namespace lawa
