namespace lawa {

template <typename T, typename Basis, typename BilinearForm>
DiagonalMatrixPreconditioner1D<T,Basis,BilinearForm>::DiagonalMatrixPreconditioner1D(const BilinearForm &_a)
    : a(_a)
{
}

template <typename T, typename Basis, typename BilinearForm>
T
DiagonalMatrixPreconditioner1D<T,Basis,BilinearForm>::operator()(XType xtype, int j, int k) const
{
    return 1./std::sqrt(fabs(a(xtype,j,k,xtype,j,k)));
}

template <typename T, typename Basis, typename BilinearForm>
T
DiagonalMatrixPreconditioner1D<T,Basis,BilinearForm>::operator()(const Index1D &index) const
{
    return DiagonalMatrixPreconditioner1D<T,Basis,BilinearForm>::operator()(index.xtype, index.j, index.k);
}

}   // namespace lawa
