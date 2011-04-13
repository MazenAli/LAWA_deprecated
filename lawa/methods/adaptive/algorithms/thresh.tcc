namespace lawa {

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index >
THRESH(const Coefficients<Lexicographical,T,Index > &v, T eta) 
{
    typedef typename Coefficients<AbsoluteValue,T,Index >::iterator it;
    Coefficients<AbsoluteValue,T,Index > temp(v.d,v.d_);
    temp = v;
    if (temp.size() > 0) {
        it lambda = temp.begin();
        T sum = 0, bound = temp.norm()*temp.norm() - eta*eta;
        do {
            sum += ((*lambda).first)*((*lambda).first);
            ++lambda;
        } while ((lambda != temp.end()) && (sum < bound));
        temp.erase(lambda, temp.end());
    }
    Coefficients<Lexicographical,T,Index > ret(v.d,v.d_);
    ret = temp;
    return ret;
}

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index >
THRESH(const Coefficients<AbsoluteValue,T,Index > &v, T eta) 
{
    typedef typename Coefficients<AbsoluteValue,T,Index >::iterator it;
    Coefficients<AbsoluteValue,T,Index > temp(v.d,v.d_);
    temp = v;
    if (temp.size() > 0) {
        it lambda = temp.begin();
        T sum = 0, bound = temp.norm()*temp.norm() - eta*eta;
        do {
            sum += ((*lambda).first)*((*lambda).first);
            ++lambda;
        } while ((lambda != temp.end()) && (sum < bound));
        temp.erase(lambda, temp.end());
    }
    Coefficients<Lexicographical,T,Index > ret(v.d,v.d_);
    ret = temp;
    return ret;
}

} // namespace lawa