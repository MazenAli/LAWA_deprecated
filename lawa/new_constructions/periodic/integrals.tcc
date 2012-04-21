#ifndef LAWA_CONSTRUCTIONS_PERIODIC_INTEGRALS_TCC
#define LAWA_CONSTRUCTIONS_PERIODIC_INTEGRALS_TCC 1

#include <lawa/math/math.h>
#include <lawa/constructions/support.h>

namespace lawa {

// reduce periodic integral to equivalent integral on R.
template <typename T>
void
_adapt_k(const Support<T> &supp1,
         const Support<T> &supp2, 
         int j1, int &k1, int j2, int &k2)
{
    if (supp1.l2>1.) {
        if (supp2.l2<supp1.l1) {
        k1 -= pow2i<T>(j1);
        }
    } else if (supp2.l2>1.) {
        if (supp1.l2<supp2.l1) {
            k2 -= pow2i<T>(j2);
        }
    } else if (supp1.l1<0.) {
        if (supp2.l1>supp1.l2) {
            k1 += pow2i<T>(j1);
        }
    } else if (supp2.l1<0) {
        if (supp1.l1>supp2.l2) {
            k2 += pow2i<T>(j2);
        }
    }
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_PERIODIC_INTEGRALS_TCC
