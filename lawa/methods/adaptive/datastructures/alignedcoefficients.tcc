namespace lawa {

template <typename T, typename Index, typename PrincipalIndex, typename AlignedIndex>
AlignedCoefficients<T,Index,PrincipalIndex,AlignedIndex>::AlignedCoefficients(void)
: map(), n1(0), n2(0)
{

}

template <typename T, typename Index, typename PrincipalIndex, typename AlignedIndex>
AlignedCoefficients<T,Index,PrincipalIndex,AlignedIndex>::AlignedCoefficients(size_t _n1, size_t _n2)
: /*map()*/ map(_n1), n1(_n1), n2(_n2)
{

}

template <typename T, typename Index, typename PrincipalIndex, typename AlignedIndex>
void
AlignedCoefficients<T,Index,PrincipalIndex,AlignedIndex>
::align_x1(const Coefficients<Lexicographical,T,Index> &coeff, short J)
{
    for (const_coeff_index_it it=coeff.begin(); it!=coeff.end(); ++it) {
        //todo: works only in 2d!! additional routine for index splitting required.
        map_prinindex_it p_prinindex=map.find((*it).first.index1);
        if (p_prinindex!=map.end()) {
            (*p_prinindex).second.operator[]((*it).first.index2) = (*it).second;
        }
        else {
            size_t tmp = std::max((size_t)pow2i<long>(J-(*it).first.index1.j+2),n2);
            Coefficients<Lexicographical,T,PrincipalIndex> coeff_x2;
            map[(*it).first.index1] = coeff_x2;
            map_prinindex_it p_prinindex=map.find((*it).first.index1);
            (*p_prinindex).second.resize(tmp);
            (*p_prinindex).second.operator[]((*it).first.index2) = (*it).second;
        }
    }
}

template <typename T, typename Index, typename PrincipalIndex, typename AlignedIndex>
void
AlignedCoefficients<T,Index,PrincipalIndex,AlignedIndex>
::unalign_x1(Coefficients<Lexicographical,T,Index> &coeff)
{
    for (const_map_prindex_it row=map.begin(); row!=map.end(); ++row) {
        for (const_coeff_aligindex_it col=(*row).second.begin(); col!=(*row).second.end(); ++col) {
            if (fabs((*col).second)>0) {
                coeff[Index2D((*row).first,(*col).first)] += (*col).second;
            }
        }
    }
}

template <typename T, typename Index, typename PrincipalIndex, typename AlignedIndex>
void
AlignedCoefficients<T,Index,PrincipalIndex,AlignedIndex>
::align_x2(const Coefficients<Lexicographical,T,Index> &coeff, short J)
{
    for (const_coeff_index_it it=coeff.begin(); it!=coeff.end(); ++it) {
        //todo: works only in 2d!! additional routine for index splitting required.
        map_prinindex_it p_prinindex=map.find((*it).first.index2);
        if (p_prinindex!=map.end()) {
            (*p_prinindex).second.operator[]((*it).first.index1) = (*it).second;
        }
        else {
            size_t tmp = std::max((size_t)pow2i<long>(J-(*it).first.index2.j+2),n2);
            Coefficients<Lexicographical,T,PrincipalIndex> coeff_x1;
            map[(*it).first.index2] = coeff_x1;
            map_prinindex_it p_prinindex=map.find((*it).first.index2);
            (*p_prinindex).second.resize(tmp);
            (*p_prinindex).second.operator[]((*it).first.index1) = (*it).second;
        }
    }
}

template <typename T, typename Index, typename PrincipalIndex, typename AlignedIndex>
std::ostream& operator<< (std::ostream &s,
                          const AlignedCoefficients<T,Index,PrincipalIndex,AlignedIndex> &alignedcoeff)
{
    s << std::endl << "AlignedCoefficients:" << std::endl;
    for (typename AlignedCoefficients<T,Index,PrincipalIndex,AlignedIndex>::const_map_prindex_it row=alignedcoeff.map.begin();
            row!=alignedcoeff.map.end(); ++row)
    {
        for (typename AlignedCoefficients<T,Index,PrincipalIndex,AlignedIndex>::const_coeff_aligindex_it col=(*row).second.begin();
                col!=(*row).second.end(); ++col)
        {
            s << (*row).first << ", " << (*col).first  << " : " << (*col).second << std::endl;
        }
    }
    return s << std::endl;
}

}   // namespace lawa
