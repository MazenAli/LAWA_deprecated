// TODO: superfluous function: replace by STL call
#ifndef LAWA_ADAPTIVE_MERGEVECTORS_H
#define LAWA_ADAPTIVE_MERGEVECTORS_H 1

#include <flens/flens.h>
#include <list>

namespace lawa {

template <typename T>
void
mergeVectors(const DenseVector<Array<T> > &in1, const DenseVector<Array<T> > &in2,
             DenseVector<Array<T> > &out)
{
    std::list<T> temp;
    for (int i=in1.firstIndex(); i<=in1.lastIndex(); ++i) {
        temp.push_back(in1(i));
    }
    for (int i=in2.firstIndex(); i<=in2.lastIndex(); ++i) {
        temp.push_back(in2(i));
    }
    temp.sort();
    temp.unique();
    out.engine().resize(temp.size(),1,0.);
    typename std::list<T>::const_iterator it=temp.begin();
    for (int i=1; it!=temp.end(); ++it, ++i) {
        out(i)=*it;
    }
}

} // namespace lawa

#endif // LAWA_ADAPTIVE_MERGEVECTORS_H
