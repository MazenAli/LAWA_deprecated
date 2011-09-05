namespace lawa {

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index >
ABSOLUTE_THRESH(const Coefficients<Lexicographical,T,Index > &v, T eta)
{
    typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
    Coefficients<Lexicographical,T,Index > ret;
    if (v.size() > 0) {
        for (const_coeff_it it=v.begin(); it!=v.end(); ++it) {
            if (fabs((*it).second) > eta) {
                ret[(*it).first] = (*it).second;
            }
        }
    }
    return ret;
}

template <typename T>
Coefficients<Lexicographical,T,Index1D>
THRESH(const Coefficients<Lexicographical,T,Index1D> &v, T eta, bool deleteBSpline, bool hp)
{
    typedef typename Coefficients<AbsoluteValue,T,Index1D >::iterator it;
    typedef typename Coefficients<Lexicographical,T,Index1D >::const_iterator const_coeff_it;

    Coefficients<AbsoluteValue,T,Index1D > temp;
    Coefficients<Lexicographical,T,Index1D > ret;
    temp = v;
    if (hp) {
        ret = THRESH(temp, eta, true);
    }
    else {
        if (temp.size() > 0) {
            it lambda = temp.begin();
            long double sum = 0L, bound = (long double)(temp.norm()*temp.norm() - eta*eta);
            do {
                sum += (long double) ((*lambda).first)*((*lambda).first);
                ++lambda;
            } while ((lambda != temp.end()) && (sum < bound));
            temp.erase(lambda, temp.end());
        }
        ret = temp;
    }

    if (!deleteBSpline) {
        for (const_coeff_it it=v.begin(); it!=v.end(); ++it) {
            if ((*it).first.xtype==XBSpline && ret.count((*it).first)==0) {
                ret[(*it).first] = (*it).second;
            }
        }
    }
    return ret;
}

template <typename T>
Coefficients<Lexicographical,T,Index2D>
THRESH(const Coefficients<Lexicographical,T,Index2D> &v, T eta, bool deleteBSpline, bool hp)
{
    typedef typename Coefficients<AbsoluteValue,T,Index2D >::iterator it;
    typedef typename Coefficients<Lexicographical,T,Index2D >::const_iterator const_coeff_it;
    Coefficients<Lexicographical,T,Index2D > ret;
    Coefficients<AbsoluteValue,T,Index2D > temp;
    temp = v;
    if (hp) {
        ret = THRESH(temp, eta, true);
    }
    else {
        if (temp.size() > 0) {
            it lambda = temp.begin();
            long double sum = 0L, bound = (long double) (temp.norm()*temp.norm() - eta*eta);
            do {
                sum += (long double) ((*lambda).first)*((*lambda).first);
                ++lambda;
            } while ((lambda != temp.end()) && (sum < bound));
            temp.erase(lambda, temp.end());
        }
        ret = temp;
    }
    if (!deleteBSpline) {
        for (const_coeff_it it=v.begin(); it!=v.end(); ++it) {
            if (    (*it).first.index1.xtype==XBSpline && (*it).first.index2.xtype==XBSpline
                 && (ret.count((*it).first)==0) ) {
                ret[(*it).first] = (*it).second;
            }
        }
    }

    return ret;
}

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index >
THRESH(const Coefficients<AbsoluteValue,T,Index > &v, T eta, bool hp)
{
    typedef typename Coefficients<AbsoluteValue,T,Index >::iterator it;
    typedef typename Coefficients<AbsoluteValue,T,Index >::reverse_iterator reverse_it;
    Coefficients<AbsoluteValue,T,Index > temp;
    temp = v;
    if (temp.size() > 0) {
        /*
         * Attention: If tol is small, then small coefficients cannot be added any more since
         * precision is not high enough.
         */
        if (!hp) {
            it lambda = temp.begin();
            long double sum = 0L, bound = (long double)(temp.norm()*temp.norm()) - (long double)eta*eta;
            do {
                sum += (long double)(((*lambda).first)*((*lambda).first));
                //std::cerr.precision(25);
                //std::cerr << (*lambda).first << " " << (*lambda).second << " : " << sum  << " " << bound << std::endl;
                ++lambda;
            } while ((lambda != temp.end()) && (sum < bound));
            temp.erase(lambda, temp.end());
        }
        else {
            reverse_it rev_lambda = temp.rbegin();
            long double sum = 0L, bound = (long double)eta*eta;
            int count=0;
            do {
                sum += (long double)(((*rev_lambda).first)*((*rev_lambda).first));
                //std::cerr.precision(25);
                //std::cerr << (*rev_lambda).first << " " << (*rev_lambda).second << " : " << sum  << " " << bound << std::endl;
                ++rev_lambda;
                ++count;
            } while ((rev_lambda != temp.rend()) && (sum <= bound));
            std::cerr << "tol = " << eta << ", count = " << count << std::endl;
            it lambda = temp.begin();
            for (int i=0; i<=(int)temp.size()-count; ++i) {
                ++lambda;
            }
            temp.erase(lambda, temp.end());
        }
    }
    Coefficients<Lexicographical,T,Index > ret;
    ret = temp;
    return ret;
}

} // namespace lawa

