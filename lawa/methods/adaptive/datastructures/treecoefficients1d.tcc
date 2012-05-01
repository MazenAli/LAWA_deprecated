namespace lawa {

template <typename T>
CoefficientsByLevel<T>::CoefficientsByLevel(void)
: j(0), map()
{
}

template <typename T>
CoefficientsByLevel<T>::CoefficientsByLevel(short _j, size_t n)
: j(_j), map() /*map(n)*/
{
    this->set(j,n);
}

template <typename T>
void
CoefficientsByLevel<T>::set(short _j, size_t n)
{
    j = _j;

    #ifdef TRONE
        if (n>4000) map.rehash(n);
    #else
        if (n>4000) map.resize(n);
    #endif

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
        map[(*it).first] += (*it).second;
    }
    return *this;
}

template <typename T>
void
CoefficientsByLevel<T>::setToZero()
{
    for (iter it=this->map.begin(); it!=this->map.end(); ++it) {
        (*it).second = 0.;
    }
}

template <typename T>
std::ostream& operator<<(std::ostream &s, const CoefficientsByLevel<T> &_coeff_by_level)
{
    s << std::endl;
    for (typename CoefficientsByLevel<T>::const_it it=_coeff_by_level.map.begin();
         it!=_coeff_by_level.map.end(); ++it) {
        s <<  (*it).first  << " : " << (*it).second << std::endl;
    }

    return s;
}



template <typename T>
TreeCoefficients1D<T>::TreeCoefficients1D(size_t n, int basis_j0)
: offset(basis_j0-1), maxTreeLevel(JMAX)
{
    for (int l=0; l<=JMAX; ++l) {
        CoefficientsByLevel<T> tmp;
        bylevel[l] = tmp;
        bylevel[l].set(l,n);
    }
}

template <typename T>
TreeCoefficients1D<T>&
TreeCoefficients1D<T>::operator=(const TreeCoefficients1D<T> &_coeff)
{
    assert(offset=_coeff.offset);
    for (int l=0; l<=maxTreeLevel; ++l) {
        this->bylevel[l] = _coeff.bylevel[l];
    }
    return *this;
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
            this->bylevel[j-1-offset].map.insert(val_type(k, (*it).second));
        }
        else {
            this->bylevel[j-offset].map.insert(val_type(k, (*it).second));
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
            this->bylevel[j-1-offset].map[k] -= (*it).second;
        }
        else {
            this->bylevel[j-offset].map[k] -= (*it).second;
        }
    }
    return *this;
}

template <typename T>
TreeCoefficients1D<T>&
TreeCoefficients1D<T>::operator-=(const TreeCoefficients1D<T> &_coeff)
{
    assert(offset=_coeff.offset);
    for (int l=0; l<=maxTreeLevel; ++l) {
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
TreeCoefficients1D<T>::operator[](short i) const
{
    // no offset here: i might be a level, then offset = 0. Or i is just the i-th stage and the
    // method calling this routine is supposed to know how work with the result vector.
    return this->bylevel[(unsigned int)i];
}

template <typename T>
CoefficientsByLevel<T>&
TreeCoefficients1D<T>::operator[](short i)
{
    return this->bylevel[(unsigned int)i];
}

template <typename T>
void
TreeCoefficients1D<T>::setToZero()
{
    for (int l=0; l<=maxTreeLevel; ++l) {
        if (this->bylevel[l].map.size()!=0) {
            for (by_level_it it=this->bylevel[l].map.begin(); it!=this->bylevel[l].map.end(); ++it) {
                (*it).second = 0.;
            }
        }
    }
}

template <typename T>
int
TreeCoefficients1D<T>::size()
{
    int ret = 0;
    for (int l=0; l<=maxTreeLevel; ++l) {
        ret += bylevel[l].map.size();
    }
    return ret;
}

template <typename T>
int
TreeCoefficients1D<T>::getMaxTreeLevel(int j0)
{
    int j=0;
    for (int l=j0-1; l<=JMAX; ++l) {
        assert(l>=0);
        if(bylevel[l].map.size()==0) break;
        j=l;
    }
    maxTreeLevel = j;
    return j;
}

template <typename T>
int
TreeCoefficients1D<T>::getMaxTreeLevel()
{
    int j=0;
    for (int l=0; l<=JMAX; ++l) {
        if(bylevel[l].map.size()==0) break;
        j=l;
    }
    maxTreeLevel = j;
    return j;
}

template <typename T>
int
TreeCoefficients1D<T>::setMaxTreeLevel(int j)
{
    maxTreeLevel = j;
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


template<typename T, typename Basis>
void
fromTreeCofficientsToCofficients(const TreeCoefficients1D<T> &tree_v,
                                 Coefficients<Lexicographical,T,Index1D> &v)
{
    for (typename CoefficientsByLevel<T>::const_it it= tree_v.bylevel[0].map.begin();
                                          it!=tree_v.bylevel[0].map.end(); ++it) {
        v[Index1D(tree_v.offset+1,(*it).first,XBSpline)] = (*it).second;
    }
    for (int i=1; i<=JMAX; ++i) {
        if (tree_v.bylevel[i].map.size()==0) break;
        for (typename CoefficientsByLevel<T>::const_it it= tree_v.bylevel[i].map.begin();
                                              it!=tree_v.bylevel[i].map.end(); ++it) {
            v[Index1D(tree_v.offset+i,(*it).first,XWavelet)] = (*it).second;
        }
    }
}

template<typename T, typename Basis>
void
fromCofficientsToTreeCofficients(const Coefficients<Lexicographical,T,Index1D> &v,
                                 TreeCoefficients1D<T> &tree_v)
{
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator    const_coeff1d_it;
    typedef typename CoefficientsByLevel<T>::val_type                           val_type;

    for (const_coeff1d_it it=v.begin(); it!=v.end(); ++it) {
        short j     = (*it).first.j;
        long  k     = (*it).first.k;
        XType xtype = (*it).first.xtype;
        if (xtype==XBSpline) {
            assert(j==tree_v.offset+1);
            tree_v.bylevel[0].map.operator[](k) =  (*it).second;
        }
        else {
            tree_v.bylevel[j-tree_v.offset].map.operator[](k) = (*it).second;
        }
    }
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
