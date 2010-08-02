namespace lawa{

template<typename T, typename Norm>
NormPreconditioner<T,Norm>::NormPreconditioner(Norm _norm)
	: norm(_norm)
{
}

template<typename T, typename Norm>
T
NormPreconditioner<T,Norm>::operator()(bool XisSpline, int j_x, int k_x,
                                       bool YisSpline, int j_y, int k_y) const
{
    return 1./sqrt(norm(XisSpline, j_x, k_x,YisSpline, j_y, k_y));
}


} // namespace lawa