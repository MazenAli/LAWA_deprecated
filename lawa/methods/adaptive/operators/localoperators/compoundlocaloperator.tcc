namespace lawa {

template <typename Index, typename FirstLocalOperator, typename SecondLocalOperator,
          typename ThirdLocalOperator,typename FourthLocalOperator>
CompoundLocalOperator<Index,FirstLocalOperator,SecondLocalOperator,ThirdLocalOperator,FourthLocalOperator>
::CompoundLocalOperator(FirstLocalOperator &_firstLocalOp, SecondLocalOperator &_secondLocalOp)
: numOfLocalOp(2),
  firstLocalOp(_firstLocalOp), secondLocalOp(_secondLocalOp),
  thirdLocalOp(_secondLocalOp), fourthLocalOp(_secondLocalOp)
{

}

template <typename Index, typename FirstLocalOperator, typename SecondLocalOperator,
          typename ThirdLocalOperator,typename FourthLocalOperator>
CompoundLocalOperator<Index,FirstLocalOperator,SecondLocalOperator,ThirdLocalOperator,FourthLocalOperator>
::CompoundLocalOperator(FirstLocalOperator &_firstLocalOp, SecondLocalOperator &_secondLocalOp,
                        ThirdLocalOperator &_thirdLocalOp)
: numOfLocalOp(3),
  firstLocalOp(_firstLocalOp), secondLocalOp(_secondLocalOp),
  thirdLocalOp(_thirdLocalOp), fourthLocalOp(_secondLocalOp)
{

}

template <typename Index, typename FirstLocalOperator, typename SecondLocalOperator,
          typename ThirdLocalOperator,typename FourthLocalOperator>
void
CompoundLocalOperator<Index,FirstLocalOperator,SecondLocalOperator,ThirdLocalOperator,FourthLocalOperator>
::eval(const Coefficients<Lexicographical,T,Index> &v, Coefficients<Lexicographical,T,Index> &Av)
{
    switch (numOfLocalOp)
    {
        case 2:
            firstLocalOp.eval( v, Av);
            secondLocalOp.eval(v, Av);
            break;
        case 3:
            firstLocalOp.eval( v, Av);
            secondLocalOp.eval(v, Av);
            thirdLocalOp.eval(v, Av);
            break;
        default:
            std::cerr << "CompoundLocalOperator not yet implemented for " << numOfLocalOp
                      << " operators. Exit." << std::endl;
            exit(1);
    }
}

template <typename Index, typename FirstLocalOperator, typename SecondLocalOperator,
          typename ThirdLocalOperator,typename FourthLocalOperator>
template <typename Preconditioner>
void
CompoundLocalOperator<Index,FirstLocalOperator,SecondLocalOperator,ThirdLocalOperator,FourthLocalOperator>
::eval(Coefficients<Lexicographical,T,Index> &v, Coefficients<Lexicographical,T,Index> &Av,
       Preconditioner &P)
{
    for (coeff_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= P[(*it).first];
    }

    switch (numOfLocalOp)
    {
        case 2:
            firstLocalOp.eval( v, Av);
            secondLocalOp.eval(v, Av);
            break;
        case 3:
            firstLocalOp.eval( v, Av);
            secondLocalOp.eval(v, Av);
            thirdLocalOp.eval(v, Av);
            break;
        default:
            std::cerr << "CompoundLocalOperator not yet implemented for " << numOfLocalOp
                      << " operators. Exit." << std::endl;
            exit(1);
    }
    for (coeff_it it=Av.begin(); it!=Av.end(); ++it) {
        (*it).second *= P[(*it).first];
    }
    for (coeff_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= 1./P[(*it).first];
    }
}

}   // namespace lawa
