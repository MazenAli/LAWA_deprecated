

namespace lawa {

template<typename T, typename Index2D, typename Basis2D>
AverageOutput2D<T,Index2D,Basis2D>::AverageOutput2D(Basis2D _basis, T xmin, T xmax, T ymin, T ymax,
		RHS<T,Index2D, SeparableRHS2D<T, Basis2D>, H1NormPreconditioner2D<T, Basis2D> > _rhs)
: basis(_basis), x_max(xmax), x_min(xmin), y_max(ymax), y_min(ymin), rhs(_rhs)
{}


template<typename T, typename Index2D, typename Basis2D>
T
AverageOutput2D<T,Index2D,Basis2D>::operator()(const CoeffVector& coeffs_u){
	T output_value;
	typename CoeffVector::const_iterator it;
	if(abs(x_max-x_min)>1e-6 && abs(y_max-y_min)>1e-6){
		T vol=(abs(x_max-x_min))*(abs(y_max-y_min));
		for(it=coeffs_u.begin(); it!=coeffs_u.end(); ++it){
			output_value+=(*it).second * rhs.operator()((*it).first);
		}
		if(vol>1e-6){
			output_value= 1/vol * output_value;
		}
		else{
			std::cout << "Fehler bei Fläche - Fehlwert 0"<< std::endl;
			output_value= 0;
		}
	}
	else{
		std::cout << "Fehler bei Abstand - Fehlwert 0"<< std::endl;
		output_value= 0;

	}
	return output_value;
}


template<typename T, typename Index2D, typename Basis2D>
T
AverageOutput2D<T,Index2D,Basis2D>::operator()(const Index2D &lambda){
	return rhs(lambda);
}

template<typename T, typename Index2D, typename Basis2D>
Coefficients<Lexicographical,T,Index2D>
AverageOutput2D<T,Index2D,Basis2D>::operator()(const IndexSet<Index2D> &Lambda){
	return rhs(Lambda);
}

template<typename T, typename Index2D, typename Basis2D>
Coefficients<Lexicographical,T,Index2D>
AverageOutput2D<T,Index2D,Basis2D>::operator()(T tol){
	return rhs(tol);
}


}		//namespace lawa
