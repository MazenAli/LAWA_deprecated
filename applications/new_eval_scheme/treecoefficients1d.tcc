namespace lawa {

template <typename T>
CoefficientsByLevel<T>&
CoefficientsByLevel<T>::operator=(const CoefficientsByLevel<T> &_coeff)
{
    this->erase(CoefficientsByLevel<T>::begin(), CoefficientsByLevel<T>::end());
    for (const_it it=_coeff.begin(); it!=_coeff.end(); ++it) {
        this->insert(val_type((*it).first,(*it).second));
    }
    return *this;
}

template <typename T>
CoefficientsByLevel<T>&
CoefficientsByLevel<T>::operator+=(const CoefficientsByLevel<T> &_coeff)
{
    for (const_it it=_coeff.begin(); it!=_coeff.end(); ++it) {
        this->operator[]((*it).first) = (*it).second;
    }
    return *this;
}

template <typename T>
std::ostream& operator<<(std::ostream &s, const CoefficientsByLevel<T> &_coeff_by_level)
{
    s << std::endl;
    for (typename CoefficientsByLevel<T>::const_iterator it=_coeff_by_level.begin();
         it!=_coeff_by_level.end(); ++it) {
        s <<  (*it).first  << " : " << (*it).second << std::endl;
    }

    return s;
}



template <typename T>
TreeCoefficients1D<T>::TreeCoefficients1D(void)
{
    for (int l=0; l<=JMAX; ++l) {
        CoefficientsByLevel<T> tmp;
        bylevel[l] = tmp;
    }
}

template <typename T>
TreeCoefficients1D<T>&
TreeCoefficients1D<T>::operator=(const Coefficients<Lexicographical,T,Index1D> &_coeff)
{
    for (const_coeff1d_it it=_coeff.begin(); it!=_coeff.end(); ++it) {
        short j     = (*it).first.j;
        long  k     = (*it).first.k;
        XType xtype = (*it).first.xtype;
        if (xtype==XBSpline) {
            this->bylevel[j-1].insert(val_type(k, (*it).second));
        }
        else {
            this->bylevel[j].insert(val_type(k, (*it).second));
        }
    }
    return *this;
}

template <typename T>
TreeCoefficients1D<T>&
TreeCoefficients1D<T>::operator-=(const Coefficients<Lexicographical,T,Index1D> &_coeff)
{
    for (const_coeff1d_it it=_coeff.begin(); it!=_coeff.end(); ++it) {
        short j     = (*it).first.j;
        long  k     = (*it).first.k;
        XType xtype = (*it).first.xtype;
        if (xtype==XBSpline) {
            this->bylevel[j-1].operator[](k) -= (*it).second;
        }
        else {
            this->bylevel[j].operator[](k) -= (*it).second;
        }
    }
    return *this;
}

template <typename T>
const CoefficientsByLevel<T>&
TreeCoefficients1D<T>::operator[](short j) const
{
    return this->bylevel[(unsigned int)j];
}

template <typename T>
CoefficientsByLevel<T>&
TreeCoefficients1D<T>::operator[](short j)
{
    return this->bylevel[(unsigned int)j];
}

template <typename T>
T
TreeCoefficients1D<T>::norm(T factor)
{
    long double norm=0.L;
    for (int l=0; l<=JMAX; ++l) {
        if (this->bylevel[l].size()!=0) {
            for (const_by_level_it it=this->bylevel[l].begin(); it!=this->bylevel[l].end(); ++it) {
                norm += std::pow((T)(*it).second,(T)factor);
            }
        }
    }
    return std::pow((T)norm,T(1)/factor);
}

template <typename T>
std::ostream& operator<<(std::ostream &s, const TreeCoefficients1D<T> &_treecoeff)
{
    for (int l=0; l<=JMAX; ++l) {
        if (_treecoeff.bylevel[l].size()!=0) {
            s << "l = " << l << " : " << std::endl;
            for (typename TreeCoefficients1D<T>::const_by_level_it it=_treecoeff.bylevel[l].begin();
                 it!=_treecoeff.bylevel[l].end(); ++it) {
                s << "  (" << l << ", " << (*it).first  << "): " << (*it).second << std::endl;
            }
        }

    }

    return s;
}

}
