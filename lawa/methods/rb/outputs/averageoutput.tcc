

namespace lawa {

template<typename T, typename Index2D, typename Basis2D>
AverageOutput2D<T,Index2D,Basis2D>::AverageOutput2D(Basis2D _basis, T xmin, T xmax, T ymin, T ymax,
		RHS<T,Index2D, SeparableRHS2D<T, Basis2D>, H1NormPreconditioner2D<T, Basis2D> > _rhs)
: basis(_basis), x_max(xmax), x_min(xmin), y_max(ymax), y_min(ymin), rhs(_rhs)
{}


template<typename T, typename Index2D, typename Basis2D>
T
AverageOutput2D<T,Index2D,Basis2D>::operator()(const Coefficients<Lexicographical, T, Index2D>& coeffs_u,
		const Coefficients<Lexicographical, T, Index2D>& basis_functions){
	T output_value;
	if(abs(x_max-x_min)>1e-6 && abs(y_max-y_min)>1e-6){
	  	std::cout << "Abstand in jeder Richtung sinnvoll definiert" <<std::endl;
		T vol=(abs(x_max-x_min))*(abs(y_max-y_min));
		for(int i=1; i<=size(coeffs_u);i++){
			output_value+=coeffs_u(i)*rhs(1,basis_functions);
		}
		if(vol>1e-6){
			std::cout << "Fläche hat positiven Wert"<< std::endl;
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

/*template<typename T, typename Index2D, typename Basis2D>
T
AverageOutput2D<T,Index2D,Basis2D>::calculate_average_output
(const Coefficients<Lexicographical, T, Index2D>& coeffs_u,
		const Coefficients<Lexicographical, T, Index2D>& basis_functions){
	T output_value;
	//hier noch Sicherheitsvorkehrungen wenn x_max-x_min<eps, ets.
	if(abs(x_max-x_min)>1e-6 && abs(y_max-y_min)>1e-6){
	  	std::cout << "Abstand in jeder Richtung sinnvoll definiert" <<std::endl;
		T vol=(abs(x_max-x_min))*(abs(y_max-y_min));
		for(int i=1; i<=size(coeffs_u);i++){
			output_value+=coeffs_u(i)*rhs(1,basis_functions);
		}
		if(vol>1e-6){ //Auslöschung? eigentlich nicht unbedingt nötig...
			std::cout << "Fläche hat positiven Wert"<< std::endl;
			output_value= 1/vol * output_value;
		}
		else{
			std::cout << "Fehler bei Fläche - Fehlwert 0"<< std::endl;
			output_value= 0; //0 als Fehlwert nehmen ?
		}
	}
	else{
		std::cout << "Fehler bei Abstand - Fehlwert 0"<< std::endl;
		output_value= 0;

	}
	return output_value;
}*/

}		//namespace lawa
