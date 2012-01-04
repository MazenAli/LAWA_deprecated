namespace lawa {

template <typename T>
CoefficientsByLevel<T>::CoefficientsByLevel(void)
: j(0), map()
{
}

template <typename T>
CoefficientsByLevel<T>::CoefficientsByLevel(short _j, size_t n)
: j(_j), map(n)
{
}

template <typename T>
void
CoefficientsByLevel<T>::set(short _j, size_t n)
{
    j = _j;
    map.resize(n);
}

template <typename T>
CoefficientsByLevel<T>&
CoefficientsByLevel<T>::operator=(const CoefficientsByLevel<T> &_coeff)
{
    map.erase(map.begin(), map.end());
    for (const_it it=_coeff.map.begin(); it!=_coeff.map.end(); ++it) {
        map.insert(val_type((*it).first,(*it).second));
    }
    return *this;
}

template <typename T>
CoefficientsByLevel<T>&
CoefficientsByLevel<T>::operator+=(const CoefficientsByLevel<T> &_coeff)
{
    for (const_it it=_coeff.map.begin(); it!=_coeff.map.end(); ++it) {
        map[(*it).first] = (*it).second;
    }
    return *this;
}

template <typename T>
std::ostream& operator<<(std::ostream &s, const CoefficientsByLevel<T> &_coeff_by_level)
{
    s << std::endl;
    for (typename CoefficientsByLevel<T>::const_iterator it=_coeff_by_level.map.begin();
         it!=_coeff_by_level.map.end(); ++it) {
        s <<  (*it).first  << " : " << (*it).second << std::endl;
    }

    return s;
}


template <typename T>
TreeCoefficients1D<T>::TreeCoefficients1D(size_t n)
{
    for (int l=0; l<=JMAX; ++l) {
        CoefficientsByLevel<T> tmp(l,n);
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
            this->bylevel[j-1].map.insert(val_type(k, (*it).second));
        }
        else {
            this->bylevel[j].map.insert(val_type(k, (*it).second));
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
            this->bylevel[j-1].map[k] -= (*it).second;
        }
        else {
            this->bylevel[j].map[k] -= (*it).second;
        }
    }
    return *this;
}

template <typename T>
TreeCoefficients1D<T>&
TreeCoefficients1D<T>::operator-=(const TreeCoefficients1D<T> &_coeff)
{
    for (int l=0; l<=JMAX; ++l) {
        if (_coeff.bylevel[l].map.size()!=0) {
            for (const_by_level_it it=_coeff.bylevel[l].map.begin(); it!=_coeff.bylevel[l].map.end(); ++it) {
                this->bylevel[l].map[(*it).first] -= (*it).second;
            }
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
        if (this->bylevel[l].map.size()!=0) {
            for (const_by_level_it it=this->bylevel[l].map.begin(); it!=this->bylevel[l].map.end(); ++it) {
                norm += std::pow((T)fabs((*it).second),(T)factor);
            }
        }
    }
    return std::pow((T)norm,T(1)/factor);
}

template <typename T>
std::ostream& operator<<(std::ostream &s, const TreeCoefficients1D<T> &_treecoeff)
{
    for (int l=0; l<=JMAX; ++l) {
        if (_treecoeff.bylevel[l].map.size()!=0) {
            s << "l = " << l << " : " << std::endl;
            for (typename TreeCoefficients1D<T>::const_by_level_it it=_treecoeff.bylevel[l].map.begin();
                 it!=_treecoeff.bylevel[l].map.end(); ++it) {
                s << "  (" << l << ", " << (*it).first  << "): " << (*it).second << std::endl;
            }
        }

    }

    return s;
}

}
