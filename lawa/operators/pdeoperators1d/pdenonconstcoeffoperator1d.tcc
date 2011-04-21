namespace lawa {

template <typename T, typename Basis>
PDENonConstCoeffOperator1D<T, Basis>::PDENonConstCoeffOperator1D(const Basis& _basis,
                                                                 Function<T>& _reaction_f,
                                                                 Function<T>& _convection_f,
                                                                 T _diffusion)
    : basis(_basis), reaction_f(_reaction_f), convection_f(_convection_f), diffusion(_diffusion),
      reaction_integral(reaction_f, basis, basis), convection_integral(convection_f, basis, basis),
      diffusion_integral(basis, basis)

{
}

template <typename T, typename Basis>
T
PDENonConstCoeffOperator1D<T, Basis>::operator()(XType xtype1, int j1, int k1,
                                              XType xtype2, int j2, int k2) const
{
    // diffusion * v_x *  u_x + convection * v * u_x + reaction * v * u
    return   diffusion * diffusion_integral(j1, k1, xtype1, 1, j2, k2, xtype2, 1)
            + convection_integral(j1, k1, xtype1, 0, j2, k2, xtype2, 1)
            + reaction_integral(j1, k1, xtype1, 0, j2, k2, xtype2, 0);
}

template <typename T, typename Basis>
T
PDENonConstCoeffOperator1D<T, Basis>::operator()(const Index1D &row_index, const Index1D &col_index)
const
{
    return this->operator()(row_index.xtype, row_index.j, row_index.k,
                                                         col_index.xtype, col_index.j, col_index.k);
}

}   // namespace lawa

