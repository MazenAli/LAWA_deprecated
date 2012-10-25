namespace lawa {

template<typename Basis, typename LocalOperator>
RHSTensorInterpolationPic2D<Basis,LocalOperator>::RHSTensorInterpolationPic2D
(const Basis &_basis, LocalOperator &_localOp, const Coefficients<Lexicographical,T,Index2D> &_u)
: basis(_basis), localOp(_localOp), u(_u)
{

}

template<typename Basis, typename LocalOperator>
void
RHSTensorInterpolationPic2D<Basis,LocalOperator>::initializeRHS
(const Coefficients<Lexicographical,T,Index2D> &_f)
{
    f = _f;
    f.setToZero();
    localOp.eval(u,f);
}

template<typename Basis, typename LocalOperator>
typename LocalOperator::T
RHSTensorInterpolationPic2D<Basis,LocalOperator>::operator()(const Index2D &index)
{
    const_coeff2d_it it = f.find(index);
    if (it!=f.end()) return (*it).second;

    std::cerr << "   RHSTensorInterpolationPic2D: extending rhs vector, " << index
              << " is missing." << std::endl;

    Coefficients<Lexicographical,T,Index2D> tmp;
    tmp = f;
    extendMultiTree(basis, tmp, f, "standard", false, true);

    f.setToZero();

    localOp.eval(u,f);

    it = f.find(index);
    if (it!=f.end()) return (*it).second;
    else {
        std::cerr << "Something went wrong in RHSTensorInterpolationPic2D<LocalOperator>" << std::endl;
        exit(1);
    }

    return 0.;
}

template<typename Basis, typename LocalOperator>
typename LocalOperator::T
RHSTensorInterpolationPic2D<Basis, LocalOperator>::operator()(XType xtype_x, int j_x, long k_x,
                                                              XType xtype_y, int j_y, long k_y)
{
    Index1D index1(j_x,k_x,xtype_x);
    Index1D index2(j_y,k_y,xtype_y);
    return this->operator()(Index2D(index1,index2));
}

}   // namespace lawa
