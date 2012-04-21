#ifndef LAWA_CONSTRUCTIONS_SUPPORT_TCC
#define LAWA_CONSTRUCTIONS_SUPPORT_TCC 1

#include <cassert>

namespace lawa {

template <typename T>
Support<T>::Support()
    : l1(0), l2(-1)
{
}

template <typename T>
Support<T>::Support(const T &a, const T &b)
    : l1(a), l2(b)
{
    assert(l1<=l2);
}

template <typename T>
Support<T>::~Support()
{
}

template <typename T>
T
Support<T>::length() const
{
    return l2-l1;
}

//------------------------------------------------------------------------------

template <typename T>
bool
inner(T x, const Support<T> &supp) {
    return (x>=supp.l1) && (x<=supp.l2);
}

template <typename T>
T
overlap(const Support<T> &supp1, const Support<T> &supp2)
{
    return std::min(supp1.l2, supp2.l2) - std::max(supp1.l1, supp2.l1);
}

template <typename T>
T
overlap(const Support<T> &supp1, const Support<T> &supp2, Support<T> &common)
{
    common.l1 = std::max(supp1.l1, supp2.l1);
    common.l2 = std::min(supp1.l2, supp2.l2);
    return common.l2 - common.l1;
}

template <typename T>
T
distance(const Support<T> &supp1, const Support<T> &supp2)
{
    return std::max(0.0, std::max(supp1.l1, supp2.l1) - std::min(supp1.l2, supp2.l2));
}

template <typename T>
T
distance(const Support<T> &supp1, const Support<T> &supp2, Support<T> &common)
{
    common.l1 = std::max(supp1.l1, supp2.l1);
    common.l2 = std::min(supp1.l2, supp2.l2);
    return std::max(0.0, common.l1 - common.l2);
}

template <typename T>
T
distance(const DenseVector<Array<T> > &singsupp1, const Support<T> &supp2)
{
    T first=singsupp1(singsupp1.firstIndex()), last=singsupp1(singsupp1.lastIndex());
    if (distance(Support<T>(first,last),supp2)>0) {
        return distance(Support<T>(first,last),supp2);
    }

    for (Integer i=singsupp1.firstIndex(); i<=singsupp1.lastIndex()-1; ++i) {
        if (singsupp1(i)<supp2.l1 && supp2.l2<singsupp1(i+1)) {
            return std::min(supp2.l1-singsupp1(i), singsupp1(i+1)-supp2.l2);
        }
        else if (supp2.l1 < singsupp1(i) && singsupp1(i) < supp2.l2) {
            return -std::min(singsupp1(i)-supp2.l1, supp2.l2-singsupp1(i));
        }
    }
    if (supp2.l1 < singsupp1(singsupp1.lastIndex()) && singsupp1(singsupp1.lastIndex()) < supp2.l2) {
        return -std::min(singsupp1(singsupp1.lastIndex())-supp2.l1, supp2.l2-singsupp1(singsupp1.lastIndex()));
    }
    return 0.0;
}

template <typename T>
T
distance(const Support<T> &supp1, const DenseVector<Array<T> > &singsupp2)
{
    return distance(singsupp2, supp1);
}

template <typename T, typename S>
Support<T>
operator+(const Support<T> &supp, S shift)
{
    return Support<T>(supp.l1+shift, supp.l2+shift);
}

template <typename S, typename T>
Support<T>
operator*(S factor, const Support<T> &supp)
{
    return Support<T>(factor*supp.l1, factor*supp.l2);
}

template <typename T>
std::ostream &
operator<<(std::ostream &out, const Support<T> &supp)
{
    out << "[" << supp.l1 << "," << supp.l2 << "]";
    return out;
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_SUPPORT_TCC
