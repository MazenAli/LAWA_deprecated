namespace lawa {

template <typename Index, typename FirstLocalOperator, typename SecondLocalOperator,
          typename ThirdLocalOperator,typename FourthLocalOperator>
CompoundLocalOperator<Index,FirstLocalOperator,SecondLocalOperator,ThirdLocalOperator,FourthLocalOperator>
::CompoundLocalOperator(const FirstLocalOperator &_firstLocalOp,
                        const SecondLocalOperator &_secondLocalOp)
: numOfLocalOp(2),
  firstLocalOp(_firstLocalOp), secondLocalOp(_secondLocalOp),
  thirdLocalOp(_secondLocalOp), fourthLocalOp(_secondLocalOp)
{

}

template <typename Index, typename FirstLocalOperator, typename SecondLocalOperator,
          typename ThirdLocalOperator,typename FourthLocalOperator>
void
CompoundLocalOperator<Index,FirstLocalOperator,SecondLocalOperator,ThirdLocalOperator,FourthLocalOperator>
::eval(const Coefficients<Lexicographical,T,Index> &v, Coefficients<Lexicographical,T,Index> &auxiliary,
       Coefficients<Lexicographical,T,Index> &Av)
{
    switch (numOfLocalOp)
    {
        case 2:
            firstLocalOp.eval( v,auxiliary,Av);
            secondLocalOp.eval(v,auxiliary,Av);
            break;
        case 3:
            firstLocalOp.eval( v,auxiliary,Av);
            secondLocalOp.eval(v,auxiliary,Av);
            break;
        default:
            std::cerr << "CompoundLocalOperator not yet implemented for " << numOfLocalOp
                      << " operators. Exit." << std::endl;
            exit(1);
    }
}

}   // namespace lawa
