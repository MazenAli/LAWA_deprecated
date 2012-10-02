namespace lawa {

template <typename Index, typename StiffnessMatrixLocalOperator, typename MassMatrixLocalOperator>
ThetaTimeStepLocalOperator<Index,StiffnessMatrixLocalOperator,MassMatrixLocalOperator>
::ThetaTimeStepLocalOperator(T _theta, T _timestep,
                             StiffnessMatrixLocalOperator &_StiffnessMatrixLocalOp)
: theta(_theta), timestep(_timestep), MassIsIdentity(true),
  StiffnessMatrixLocalOp(_StiffnessMatrixLocalOp), MassMatrixLocalOp(_StiffnessMatrixLocalOp)
{
    if (theta<0. || theta>1.) {
        std::cerr << "Non-admissible theta value: " << theta << std::endl;
        exit(1);
    }
    if (timestep < 0.) {
        std::cerr << "Non-admissible timestep value: " << timestep << std::endl;
        exit(1);
    }
}

template <typename Index, typename StiffnessMatrixLocalOperator, typename MassMatrixLocalOperator>
ThetaTimeStepLocalOperator<Index,StiffnessMatrixLocalOperator,MassMatrixLocalOperator>
::ThetaTimeStepLocalOperator(T _theta, T _timestep,
                             StiffnessMatrixLocalOperator &_StiffnessMatrixLocalOp,
                             MassMatrixLocalOperator &_MassMatrixLocalOp)
: theta(_theta), timestep(_timestep), MassIsIdentity(false),
  StiffnessMatrixLocalOp(_StiffnessMatrixLocalOp), MassMatrixLocalOp(_MassMatrixLocalOp)
{
    if (theta<0. || theta>1.) {
        std::cerr << "Non-admissible theta value: " << theta << std::endl;
        exit(1);
    }
    if (timestep < 0.) {
        std::cerr << "Non-admissible timestep value: " << timestep << std::endl;
        exit(1);
    }
}

template <typename Index, typename StiffnessMatrixLocalOperator, typename MassMatrixLocalOperator>
void
ThetaTimeStepLocalOperator<Index,StiffnessMatrixLocalOperator,MassMatrixLocalOperator>
::setThetaTimeStepParameters(T _theta, T _timestep)
{
    if (theta<0. || theta>1.) {
        std::cerr << "Non-admissible theta value: " << theta << std::endl;
        exit(1);
    }
    if (timestep < 0.) {
        std::cerr << "Non-admissible timestep value: " << timestep << std::endl;
        exit(1);
    }
    theta = _theta;
    timestep = _timestep;
}

template <typename Index, typename StiffnessMatrixLocalOperator, typename MassMatrixLocalOperator>
void
ThetaTimeStepLocalOperator<Index,StiffnessMatrixLocalOperator,MassMatrixLocalOperator>
::evalM(Coefficients<Lexicographical,T,Index> &v, Coefficients<Lexicographical,T,Index> &Mv,
        const char* evalType)
{
    if (MassIsIdentity) Mv += v;
    else                MassMatrixLocalOp.eval(v, Mv, evalType);
    return;
}

template <typename Index, typename StiffnessMatrixLocalOperator, typename MassMatrixLocalOperator>
void
ThetaTimeStepLocalOperator<Index,StiffnessMatrixLocalOperator,MassMatrixLocalOperator>
::evalA(Coefficients<Lexicographical,T,Index> &v, Coefficients<Lexicographical,T,Index> &Av,
        const char* evalType)
{
    StiffnessMatrixLocalOp.eval(v, Av, evalType);
    return;
}

template <typename Index, typename StiffnessMatrixLocalOperator, typename MassMatrixLocalOperator>
template <typename Preconditioner>
void
ThetaTimeStepLocalOperator<Index,StiffnessMatrixLocalOperator,MassMatrixLocalOperator>
::eval(Coefficients<Lexicographical,T,Index> &v, Coefficients<Lexicographical,T,Index> &Av,
       Preconditioner &P, const char* evalType)
{
    for (coeff_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= P((*it).first);
    }

    StiffnessMatrixLocalOp.eval(v, Av, evalType);
    Av *= (theta*timestep);
    if (MassIsIdentity) Av += v;
    else                MassMatrixLocalOp.eval(v, Av, evalType);

    for (coeff_it it=Av.begin(); it!=Av.end(); ++it) {
        (*it).second *= P((*it).first);
    }
    for (coeff_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= 1./P((*it).first);
    }
}


}   // namespace lawa
