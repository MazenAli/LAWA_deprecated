#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_MRA_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_PRIMAL_MRA_TCC 1

// Since the primal MRA for the interval construction of Dijkema is
// identical to the one of the Primbs construction, we just reuse that one.
// Note: this is not a very nice solution. It depends on the fact, that the
//       Primbs MRA implementation does NOT include other files where this
//       macro substitution could be really dangerous!
#define Primbs Dijkema
#include <lawa/constructions/interval/primbs/primal/mra.tcc>
#undef Primbs
#undef LAWA_CONSTRUCTIONS_INTERVAL_PRIMBS_PRIMAL_MRA_TCC

#endif //LAWA_INTERVAL_DIJKEMA_PRIMAL_MRA_TCC
