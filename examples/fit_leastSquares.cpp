#include <iostream>

#include <LINALG/matrix.hpp>
#include <LINALG/solve.hpp>

//int d, const std::string  filePath, int n
int main (int argc, char* argv[])
{
	int d = stoi(std::string(argv[1])); //partition size
	dMatrix data = getDataSet(argv[2]);	//data set
	int n = stoi(std::string(argv[3])); //polinomial degree

	double minx = std::numeric_limits<double>::infinity();
	double maxx = -1*std::numeric_limits<double>::infinity();
	for (int i = 0; i < data.getRows(); i++)
	{
		minx = (data[i][0] < minx) ? data[i][0] : minx;
		maxx = (data[i][0] > maxx) ? data[i][0] : maxx;
	}

	// minx = -15;
	// maxx = 15;

	dMatrix A = generatePhi(n,data);
	dVector b = data.col(1);

	dMatrix coefficients = fitLeastSquares(A,b);
	
	if (coefficients.getType() == dVector::MatrixType::empty)
	{
		std::cout<<"Coefficients smaller than tolerance"<<std::endl;
		std::cout<<"Not posible to fit!"<<std::endl;
		return 0;
	}

	std::cout<<"Coefficient vector: "<<coefficients.T()<<"Error: "<<MSE(b,A*coefficients)<<std::endl;
	double* partition = (double*)alloca(d * sizeof(double));
	for (int i = 0; i<d; i++)
		partition[i] = minx + (maxx-minx)*i/(d);

	exportPolinomial(n, coefficients.T()[0], d, partition, "test_out.txt");
	std::cout<<"Polinomial (x,y) exported successfully!"<<std::endl;
	return 0;
}

