#include <iostream>
#include "MATRIX/matrix.hpp"
#include "DIFF/diff_solve.hpp"

int main()
{
	auto f = [](dVector v)->double {
		double x = v[0][0], y = v[1][0];
		// inserte función de prueba

		return 0;
	};

	plot3d("surface.txt",f,-5,5,-5,5);

	auto df = [](dVector v)->dVector {
		double x = v[0][0], y = v[1][0];
		dVector g(2,1);
		//inserte derivada de la función de prueba

		return g;
	};

	dVector v(2,1);
	v[0][0] = 2;
	v[1][0] = 1;

	double epsilon = std::sqrt(getMachineEpsilon<double>());

	std::cout<<"f(x,y) = x^2 + y^2"<<std::endl;
	std::cout<<"step size = 0.1"<<std::endl;
	std::cout<<"======================"<<std::endl;
	{
		dVector autoargmin = autograd(f,v,0.1,epsilon,false);
		dVector argmin = optimize(df,v,0.1,epsilon, false);
		std::cout<<"autoargmin = "<<autoargmin.T()<<std::endl;
		std::cout<<"argmin = "<<argmin.T()<<std::endl;
	}
	std::cout<<"step size = 0.5"<<std::endl;
	std::cout<<"======================"<<std::endl;
	{
		dVector autoargmin = autograd(f,v,0.5,epsilon,false);
		dVector argmin = optimize(df,v,0.5,epsilon,false);
		std::cout<<"autoargmin = "<<autoargmin.T()<<std::endl;
		std::cout<<"argmin = "<<argmin.T()<<std::endl;
	}

	return 0;
}
