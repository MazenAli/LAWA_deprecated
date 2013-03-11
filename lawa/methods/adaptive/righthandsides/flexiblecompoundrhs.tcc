namespace lawa {

template <typename T, typename Index, typename RHSType>
FlexibleCompoundRhs<T, Index,RHSType>::FlexibleCompoundRhs(std::vector<RHSType*>& _rhsvec)
 : rhsvec(_rhsvec), active_comp(rhsvec.size())
{
	for(int i=0; i < active_comp.size(); ++i){
		active_comp[i]=i;
	}
}

template <typename T, typename Index, typename RHSType>
T
FlexibleCompoundRhs<T, Index,RHSType>::operator()(const Index &index)
{
    T val = 0.;
//    for(auto& r : rhsvec){
//    	val += (*r)(index);
//    }
    for(int i : active_comp){
    	val += (*rhsvec[i])(index);
    }

    return val;
}

template <typename T, typename Index, typename RHSType>
void
FlexibleCompoundRhs<T, Index,RHSType>::set_active_comp(int i)
{
	if(i < 0){
		active_comp.resize(rhsvec.size());
		for(int i=0; i < active_comp.size(); ++i){
			active_comp[i]=i;
		}
	}
	else{
		assert(i < rhs_vec.size());
		active_comp.resize(1);
		active_comp[0] = i;
	}
}

}   // namespace lawa

