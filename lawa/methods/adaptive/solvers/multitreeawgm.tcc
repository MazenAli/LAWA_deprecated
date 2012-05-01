namespace lawa {

template <typename Index, typename LocalOperator, typename RHS, typename Preconditioner>
MultiTreeAWGM<Index,LocalOperator,RHS,Preconditioner>::
MultiTreeAWGM(LocalOperator &_A, RHS &_F, Preconditioner &_P)
: A(_A), F(_F), P(_P), alpha(0.5), gamma(0.1)
{

}

template <typename Index, typename LocalOperator, typename RHS, typename Preconditioner>
void
MultiTreeAWGM<Index,LocalOperator,RHS,Preconditioner>::setParameters(T _alpha, T _gamma)
{
    alpha = _alpha;
    gamma = _gamma;
}

template <typename Index, typename LocalOperator, typename RHS, typename Preconditioner>
void
MultiTreeAWGM<Index,LocalOperator,RHS,Preconditioner>::
solve(Coefficients<Lexicographical,T,Index> &u, T _eps, const char *filename,
      int NumOfIterations, T EnergyNorm)
{
    Coefficients<Lexicographical,T,Index2D> r(SIZEHASHINDEX2D), p(SIZEHASHINDEX2D),
                                                Ap(SIZEHASHINDEX2D);
}

}   // namespace lawa
