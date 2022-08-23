#include "intg_solve.hpp"

double RTrapezoidal(std::function<double(double)>f, double a, double b, int n, double* storage)
{
	//No voy a necesitar esta función. Hacer esto solo sería mas lento.
	return 0;
}

dMatrix getRombergMat(std::function<double(double)> f, double a, double b, int n)
{
	dMatrix rStore(n+1,n+1,dMatrix::MatrixType::zeros);

	rStore[0][0] = ((b-a)/2) * (f(a)+f(b));
	
	// A bottom up aproach is better;
	for (int i = 1; i <= n; i++)
	{
		double h = (b-a)/pow(2.0,i);
		double sum = 0;
		double sumLimit = pow(2,i-1);
		for (int k=1; k <= sumLimit; k++) sum+=f(a + (2*k-1)*h);
		rStore[i][0] = 0.5*rStore[i-1][0]+h*sum;
	}

	for (int j = 1; j <= n; j++)
	{
		for (int i = j; i <= n; i++)
		{
			rStore[i][j] = rStore[i][j-1] + 1/(pow(4.0,j)-1) * (rStore[i][j-1]-rStore[i-1][j-1]);
		}
	}

	return rStore;
}