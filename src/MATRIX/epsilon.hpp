// Debug only
//------------------------
// #ifndef iostream
// #include <iostream>
// #endif
//------------------------

#define LEPSILON epsilon

//epislonType function out of the PDF
template <class dT>
dT getMachineEpsilon()
{
	//Using static cast on every constant to esnure correct type of operation
	dT epsilon = static_cast<dT>(1);
	dT unidad = static_cast<dT>(1);

	dT valor = unidad + epsilon;

	while (valor>unidad)
	{
		epsilon/=static_cast<dT>(2);
		valor = unidad + epsilon;
	}

	epsilon*=static_cast<dT>(2);

	return epsilon;
}
