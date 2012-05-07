namespace lawa {

template <typename Index, typename FirstLocalOperator, typename SecondLocalOperator,
          typename ThirdLocalOperator,typename FourthLocalOperator>
CompoundLocalOperator<Index,FirstLocalOperator,SecondLocalOperator,ThirdLocalOperator,FourthLocalOperator>
::CompoundLocalOperator(FirstLocalOperator &_firstLocalOp, SecondLocalOperator &_secondLocalOp,
                        size_t hashMapSize)
: numOfLocalOp(2),
  firstLocalOp(_firstLocalOp), secondLocalOp(_secondLocalOp),
  thirdLocalOp(_secondLocalOp), fourthLocalOp(_secondLocalOp),
  auxiliary(hashMapSize)
{

}

template <typename Index, typename FirstLocalOperator, typename SecondLocalOperator,
          typename ThirdLocalOperator,typename FourthLocalOperator>
void
CompoundLocalOperator<Index,FirstLocalOperator,SecondLocalOperator,ThirdLocalOperator,FourthLocalOperator>
::eval(const Coefficients<Lexicographical,T,Index> &v, Coefficients<Lexicographical,T,Index> &Av)
{
    auxiliary = Av;
    auxiliary.setToZero();

    switch (numOfLocalOp)
    {
        case 2:
            firstLocalOp.eval( v,auxiliary,Av);
            //std::cerr << "Av = " << Av << std::endl;
            auxiliary.setToZero();
            secondLocalOp.eval(v,auxiliary,Av);
            //std::cerr << "Av = " << Av << std::endl;
            break;
        case 3:
            firstLocalOp.eval( v,auxiliary,Av);
            auxiliary.setToZero();
            secondLocalOp.eval(v,auxiliary,Av);
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
    auxiliary = Av;
    auxiliary.setToZero();

    for (coeff_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= P[(*it).first];
    }

    switch (numOfLocalOp)
    {
        case 2:
            firstLocalOp.eval( v,auxiliary,Av);
            //std::cerr << "Av = " << Av << std::endl;
            auxiliary.setToZero();
            secondLocalOp.eval(v,auxiliary,Av);
            //std::cerr << "Av = " << Av << std::endl;
            break;
        case 3:
            firstLocalOp.eval( v,auxiliary,Av);
            auxiliary.setToZero();
            secondLocalOp.eval(v,auxiliary,Av);
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
