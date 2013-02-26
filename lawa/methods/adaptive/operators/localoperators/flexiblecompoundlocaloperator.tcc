namespace lawa {

template <typename Index, typename LocalOperatorType>
FlexibleCompoundLocalOperator<Index, LocalOperatorType>::
FlexibleCompoundLocalOperator(std::vector<LocalOperatorType*>& _localops) :
	localops(_localops) {}

template <typename Index, typename LocalOperatorType>
void
FlexibleCompoundLocalOperator<Index, LocalOperatorType>::
eval(const Coefficients<Lexicographical,T,Index> &v,
	 Coefficients<Lexicographical,T,Index> &Av)
{
	for(auto &op : localops){
		op->eval(v, Av);
	}
}

template <typename Index, typename LocalOperatorType>
template <typename Preconditioner>
void
FlexibleCompoundLocalOperator<Index, LocalOperatorType>::
eval(Coefficients<Lexicographical,T,Index> &v,
	 Coefficients<Lexicographical,T,Index> &Av, Preconditioner &P)
{
	for(auto& el : v){
		el.second *= P(el.first);
	}

	for(auto& op : localops){
		op->eval(v, Av);
	}

	for(auto& el : Av){
		el.second *= P(el.first);
	}
	for(auto& el : v){
		el.second /= P(el.first);
	}
}

template <typename Index, typename LocalOperatorType>
template <typename RightPrec, typename LeftPrec>
void
FlexibleCompoundLocalOperator<Index, LocalOperatorType>::
eval(Coefficients<Lexicographical,T,Index> &v,
	 Coefficients<Lexicographical,T,Index> &Av, RightPrec& rightP, LeftPrec& leftP)
{
	for(auto& el : v){
		el.second *= rightP(el.first);
	}

	for(auto& op : localops){
		op->eval(v, Av);
	}

	for(auto& el : Av){
		el.second *= leftP(el.first);
	}
	for(auto& el : v){
		el.second /= rightP(el.first);
	}
}


template <typename Index, typename LocalOperatorType>
void
FlexibleCompoundLocalOperator<Index, LocalOperatorType>::
eval(Coefficients<Lexicographical,T,Index> &v,
	 Coefficients<Lexicographical,T,Index> &Av,
	 Coefficients<Lexicographical,T,Index>& rightP,
	 Coefficients<Lexicographical,T,Index>& leftP)
{
	for(auto& el : v){
		el.second *= rightP[el.first];
	}

	for(auto& op : localops){
		op->eval(v, Av);
	}

	for(auto& el : Av){
		el.second *= leftP[el.first];
	}
	for(auto& el : v){
		el.second /= rightP[el.first];
	}
}

} // namespace lawa

