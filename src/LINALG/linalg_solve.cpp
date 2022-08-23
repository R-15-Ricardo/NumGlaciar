#include "linalg_solve.hpp"

// - - - - - - - - Substitutions - - - - - - - - 

dVector forwardSubstitution(const dMatrix& L, const dVector& b, double tolerance)
{
	// LOG("forwardSubstitution()");
	int n = b.getRows();
	dVector x(n,1,dMatrix::MatrixType::colvector);

	for (int i = 0; i < n; i++)
	{
		double sum = 0;
		for (int j = 0; j < i; j++) sum += L[i][j]*x[j][0];

		if (std::abs(L[i][i])<tolerance)
			return dVector();

		x[i][0] = (b[i][0]-sum)/(L[i][i]);
	}

	return x;
}

dVector backwardSubstitution(const dMatrix& U, const dVector& b, double tolerance)
{
	// LOG("backwardSubstitution()");
	int n = b.getRows();
	dVector x(n,1,dMatrix::MatrixType::colvector);

	for (int i = n-1; i >= 0; i--)
	{
		double sum = 0;
		for (int j = i+1; j < n; j++) sum += U[i][j]*x[j][0];

		if (std::abs(U[i][i])<tolerance)
			return dVector();
			
		x[i][0] = (b[i][0]-sum)/(U[i][i]);
	}

	return x;
}

// - - - - - - - - Solvers - - - - - - - - 
// - - - - - - - - - - - explicit - - - - - - - - - - -
dVector genSolLU (const dMatrix& L, const dMatrix& U, const dVector& b, unsigned int* swapWith, double tolerance)
{	
	// LOG("genSolLU()");
	unsigned int n = L.getRows();
	
	dVector b_hat(n,1);

	for (int i = 0; i < n; i++)
		b_hat[i][0] = b[swapWith[i]][0];

	dVector y = forwardSubstitution(L,b_hat,tolerance);
	if (y.getType() == dVector::MatrixType::empty)
		return dVector();

	dVector x = backwardSubstitution(U,y,tolerance);
	if (x.getType() == dVector::MatrixType::empty)
		return dVector();
	
	return x;
}

dVector solFactLU (const dMatrix A, const dVector b)
{
	// LOG("solFactLU()");
	double tolerance = pow(getMachineEpsilon<double>(),2.0/3.0);
	unsigned int n = A.getRows();
	dMatrix L(n,n),U(n,n);

	unsigned int p[n];

	if (factLU(A,L,U,p,tolerance))
		//Para el ayudante: Los objetos son destruidos al salir del scope. ~dVector (a.k.a. ~Matrix<double> libera memoria de manera normal) 
		return dVector();

	dVector x = genSolLU(L,U,b,p,tolerance);
	if (x.getType() == dVector::MatrixType::empty)
		return dVector();

	return x;
}

dVector solFactChol (const dMatrix A, const dVector b, double tolerance)
{
	// LOG("solFactChol()");
	unsigned int n = A.getRows();

	dMatrix L = factChol(A,tolerance);
	if (L.getType() == dMatrix::MatrixType::empty)
		return dVector();
	
	unsigned int p[n];
	for (int i = 0; i < n; i++) p[i] = i;

	dVector x = genSolLU(L,L.T(),b,p,tolerance);
	if (x.getType() == dMatrix::MatrixType::empty)
		return dVector();
	
	return x;
}

dVector solTriD_Thomas(const dMatrix& B, const dVector& d, double tolerance)
{
	int n = B.getRows();

	dMatrix x(n,1);

	dMatrix bar(n,3);

	// b bar
	bar[0][0] = B[0][1];
	// c bar
	bar[0][1] = B[0][2];
	// d bar
	bar[0][2] = d[0][0];

	for (int i = 1; i < n; i++)
	{ 
		// b bar
		bar[i][0] = bar[i-1][0]*B[i][1] - B[i][0]*bar[i-1][1];
		// c bar
		bar[i][1] = bar[i-1][0]*B[i][2];
		// d bar
		bar[i][2] = bar[i-1][0]*d[i][0] - B[i][0]*bar[i-1][2];

		//LOG(bar[i][0]<<" "<<bar[i][1]<<" "<<bar[i][2]);
	}

	if (std::abs(bar[n-1][0]) < tolerance)
		return dVector();

	x[n-1][0] = bar[n-1][2]/bar[n-1][0];

	for (int i = n-2; i>=0; i--)
	{
		if (std::abs(bar[i][0]) < tolerance)
			return dVector();
		
		x[i][0] = (bar[i][2]-bar[i][1]*x[i+1][0])/bar[i][0];
	}

	return x;
}

dVector fitLeastSquares (const dMatrix& A, const dVector& b)
{
	dMatrix FitMatrix = A.T()*A;
	dVector fitVector = A.T()*b;

	double tolerance = pow(getMachineEpsilon<double>(), 2.0/3.0);

	dVector fit = solFactChol(FitMatrix, fitVector, tolerance);
	
	// No-solution exception returned by solFactChol(): MatrixType::empty
	return fit;
}

// - - - - - - - - - - - iterative - - - - - - - - - - -

int solTriD_GS (const dMatrix& B, const dVector& b, dVector& x, unsigned int maxIter, double tolerance, double& error)
{
	int n = B.getRows();
	dVector errorVector(n,1);

	int reportInterval = maxIter/10;

	int iter;
	for (iter = 1; iter <= maxIter; iter++)
	{
		//Update vector
		for (int i = 0; i<n; i++)
		{
			double t1 = 0, t2 = 0;
			if (std::abs(B[i][1]) < tolerance)
			{
				error = std::numeric_limits<double>::infinity();
				return -1;
			}
	
			t1 = (0 <= i-1) ? B[i][0]*x[i-1][0] : 0;
			t2 = (i+1 < n) ? B[i][2]*x[i+1][0] : 0;

			x[i][0] = (b[i][0] - t1 - t2)/(B[i][1]);
		}
		//Calculate error
		errorVector.setType(dVector::MatrixType::zeros);
		for (int i = 0; i<n; i++)
		{
			for (int j = 0; j < 3; j++)
			{	
				double aMatFactor = (0<=i+j-1 && i+j-1<n) ? B[i][j]*x[i+j-1][0] : 0;
				errorVector[i][0] += aMatFactor;
			}

			errorVector[i][0] = errorVector[i][0] - b[i][0];
		}

		error = normL2(errorVector);

		if ((iter-1)%reportInterval == 0)
		{
			std::cout<<"solveTriD():[ Iteration: "<<iter<<" | error: "<<error<<" ]"<<std::endl;
			std::cout<<"x("<<iter<<") = "<<x.T()<<std::endl;
		}

		if (error < tolerance)
			return iter;
	}
	return iter-1;
}

dVector approxSolSVD(const dVector& b, const dMatrix& U, const dMatrix& V, std::vector<double>& s, unsigned int k)
{
	dVector x(U.getRows(), 1, dVector::MatrixType::zeros);

	for (int i = 0; i < k; i++)
	{
		x = x + (U.col(i).T()*b)[0][0]/s[i] * V.col(i);
	}

	return x;
}


// - - - - - - - - Factorizations - - - - - - - - 
bool factLU(const dMatrix A,dMatrix& L,dMatrix& U, unsigned int* p, double tolerance)
{
	// LOG("factLU()");
	unsigned int n = A.getRows();
	L.setType(dMatrix::MatrixType::identity);
	U = A;
	for (int i = 0; i < n; i++)  p[i] = i;

	//Indecies as per PDF
	for (int k = 0; k < n; k++)
	{
		unsigned int rIndex=0;
		double maxVal=0;
		for (int i = k; i < n; i++)
		{
			if (std::abs(U[i][k]) > maxVal) 
			{
				maxVal = U[i][k];
				rIndex = i;
			}
		}

		if (std::abs(U[rIndex][k]) < tolerance)
			return 1;
		
		if (rIndex != k)
		{
			U.swapRows(k,rIndex);

			int aux = p[k];
			p[k] = p[rIndex];
			p[rIndex] = aux;
		
			if (k>0) 
			{
				for (int j = 0; j < k; j++)
				{
					double aux_1 = L[k][j];
					L[k][j] = L[rIndex][j]; 
					L[rIndex][j] = aux_1;
				}
			}
		}

		for (int i = k+1; i < n; i++)
		{
			L[i][k] = U[i][k]/U[k][k];
			for (int j = k; j < n;j++) 
			{
				U[i][j] = U[i][j] - L[i][k] * U[k][j];
				if (std::abs(U[i][j]) < tolerance) U[i][j] = 0;
			}
		}
	}
	return 0;
}

dMatrix factChol(const dMatrix A, double tolerance)
{
	// LOG("factChol()");
	unsigned int n = A.getRows();
	dMatrix L(n,n,dMatrix::MatrixType::zeros);

	for (int j = 0; j < n; j++)
	{
		double sum = 0;
		for (int k = 0; k < j; k++) sum += L[j][k]*L[j][k];
		double rad = A[j][j] - sum;

		if (rad < 0) return dMatrix();

		L[j][j] = sqrt(rad);

		for (int i = j+1; i < n; i++)
		{
			sum = 0;
			for (int k = 0; k < j; k++) sum += L[i][k]*L[j][k];

			if (L[j][j] < tolerance) return dMatrix();
			
			L[i][j] = (A[i][j] - sum)/L[j][j];
		}
	}

	return L;
}

std::vector<double> iterativeSVD (dMatrix A, double tolerance, unsigned int maxIterations, dMatrix& U, dMatrix& V)
{
	double MachineEpsilon = getMachineEpsilon<double>();
	unsigned int k = 0, F = 1;
	unsigned int m = A.getRows(), n = A.getCols();
	
	std::vector< std::pair<double, int> > singularValues;

	V = dMatrix(n,n, dMatrix::MatrixType::identity);
	U = dMatrix(m,n, dMatrix::MatrixType::identity);

	while (k < maxIterations && F > 0)
	{
		k++;
		F = 0;

		int logMaster = (static_cast<int>((n-1)/10) > 0) ? static_cast<int>((n-1)/10) : 1;
		std::cout<<"Iteracion "<<k<<"/"<<maxIterations<<std::flush;

		for (int i = 0; i < n-1; i++)
		{

			if (i % logMaster == 0 ) std::cout<<"."<<std::flush;

			for (int j = i+1; j < n; j++)
			{
				double alpha = (A.col(i).T()*A.col(i))[0][0];
				double gamma = (A.col(j).T()*A.col(j))[0][0];
				double beta = (A.col(i).T()*A.col(j))[0][0];

				
				if (alpha * gamma >  MachineEpsilon && std::abs(beta) > tolerance*alpha*gamma)
				{
					F = 1;

					double c, s; 
					if (beta != 0)
					{
						double eta = (gamma-alpha)/(2*beta);
						double t = 1/(std::abs(eta)+sqrt(1+eta*eta));
						if (eta < 0) t *= -1;

						c = 1 / sqrt(1+t*t);
						s = t*c;
					}
					else 
					{
						c = 1;
						s = 0;
					}

					dVector a,b;

					a = A.col(i);
					b = A.col(j);

					A.insertCol(c*a - s*b, i);
					A.insertCol(s*a + c*b, j);

					a = V.col(i);
					b = V.col(j);


					V.insertCol(c*a - s*b, i);
					V.insertCol(s*a + c*b, j);

				}
			}
		}
		std::cout<<std::endl;
	}

	for (int i = 0; i < n; i++)
	{
		double s; 
		s = normL2(A.col(i));
		singularValues.push_back(std::pair<double, int>(s,i));
	}

	std::sort(singularValues.begin(), singularValues.end());
	std::reverse(singularValues.begin(), singularValues.end());

	dMatrix auxA(m,n, dMatrix::MatrixType::zeros);
	dMatrix auxV(n,n, dMatrix::MatrixType::zeros);

	for (int i = 0; i < n; i++)
	{
		std::pair<double, int> p = singularValues[i];
		auxA.insertCol(A.col(p.second),i);
		auxV.insertCol(V.col(p.second),i);
	}

	A = auxA;
	V = auxV;

	for (int i = 0; i < n; i++)
	{
		dVector u;
		u = (1/singularValues[i].first) * A.col(i);
		U.insertCol(u,i);
	}

	std::vector<double> sv;
	for (auto p : singularValues)
		sv.push_back(p.first);


	return sv;
}

// - - - - - - - - Helper functions - - - - - - - - 

dMatrix reduceTriDMatrix(const dMatrix& A)
{
	int n = A.getRows();

	dMatrix B(n,3);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < 3; j++)
			B[i][j] = (i+j-1 < 0 || i+j-1 >= n) ? 0 : A[i][i+j-1];

	return B;
}

dMatrix getDataSet (const char filePath[])
{
	using namespace std;

	fstream f;
	string fileLine;
	f.open(filePath, ios::in);

	int dataCount = 0;
	if (f.is_open())
	{
		while(getline(f,fileLine))
			dataCount++;
	}
	f.close();

	f.open(filePath, ios::in);
	dMatrix data(dataCount,2);

	fileLine = "";
	int i = 0;	
	if (f.is_open())
	{
		while(getline(f,fileLine))
		{
			fileLine+=' ';
			string num = ""; 
			bool xory = 0;
			for (char c: fileLine) 
			{ 
				if (c == ' ') 
				{
					if (xory == 0) data[i][0] = stod(num);
					else data[i][1] = stod(num);

					xory = 1;
					num = "";
				}
				else
					num += c;
			}
			i++;
		}
	}
	f.close();

	return data;
}

dMatrix generatePhi (int degree, const dMatrix& data)
{
	unsigned int pointsCount = data.getRows();
	
	dMatrix phi(pointsCount, degree+1, dMatrix::MatrixType::identity);

	for (int i = 0; i < pointsCount; i++)
		for (int j = degree; j >= 0; j--)
			phi[i][degree-j] = pow(data[i][0],j);

	return phi;
}

double evaluateHorner(int degree, double* coefficients, double x)
{
	double* b = (double*)alloca(degree * sizeof(double));
	b[0] = coefficients[0];
	for (int i = 1; i<degree; i++)
		b[i] = coefficients[i]+b[i-1]*x;

	return coefficients[degree] + b[degree-1]*x;
}

double evaluateSpline (double x, const dMatrix& data, const std::vector<double> M)
{
	unsigned int n = data.getRows()-1;
	for (int i = 1; i <= n; i++)
	{
		if (data[i-1][0] <= x && x <= data[i][0])
		{	
			double h = data[i][0] - data[i-1][0];
			double C_hat = data[i-1][1] - M[i-1] * pow(h,2)/6;
			double C = (data[i][1] - data[i-1][1])/h - h/6 * (M[i] - M[i-1]);

			double spline = M[i-1] * (pow(data[i][0] - x,3))/(6*h) + M[i] * (pow(x-data[i-1][0],3)/(6*h)) + C*(x-data[i-1][0]) + C_hat;
			return spline;
		}
	}

	return 0;
}

void exportPolinomial(int degree, double* coefficients, int partitionSize, double* partition, const std::string outPath)
{
	using namespace std;
	fstream f;
	f.open(outPath.c_str(),ios::out);

	if (f.is_open())
	{
		for (int i = 0; i<partitionSize; i++)
			f<<partition[i]<<" "<<evaluateHorner(degree, coefficients, partition[i])<<endl;
	}
	f.close();
}

void exportSpline (double a, double b, const dMatrix& data, const std::vector<double> M, const std::string& outPath)
{
	using namespace std;
	fstream file;
	file.open(outPath.c_str(),ios::out);

	if (file.is_open())
	{
		for (int i = 0; i<=100; i++)
		{
			double z = a +  i*(b-a)/100;
			double f = z + (1+z)*sin(z*z) + 2*cos(6*z);
			double s = evaluateSpline(z,data,M);

			file<<z<<" "<<f<<" "<<s<<std::endl;
		}
	}
	std::cout<<"Spline exported successfully!"<<std::endl;
	file.close();
}

double matrixConditionNumber(const dMatrix& A)
{
	std::vector<double> sv;
	dMatrix U,V;
	sv = iterativeSVD(A,sqrt(getMachineEpsilon<double>()),100,U,V);
	
	return sv.front()/sv.back();
}

std::vector<double> getD2CSplineCoef(const dMatrix& data)
{
	unsigned int n = data.getRows() - 1;

	dMatrix B(n-1, 3);
	dVector d(n-1, 1);

	//data[i][0] = x_i
	//data[i][1] = y_i = f_i 

	// calculate d's
	for (int i = 1; i <= n-1; i++)
	{
		double hi  = data[i][0] - data[i-1][0];
		double hi1 = data[i+1][0] - data[i][0];
		d[i-1][0] = 6.0/(hi+hi1)*((data[i+1][1]-data[i][1])/hi1 - (data[i][1]-data[i-1][1])/hi);
	}

	// calculate b
	for (int i = 1; i <= n-1; i++)
	{
		double h0, h1; // h0 = h_i | h1 = h_{i+1}
		h0 = data[i][0] - data[i-1][0];
		h1 = data[i+1][0] - data[i][0];

		double mu, lambda;

		mu = h0/(h0+h1);
		lambda = h1/(h0+h1);

		B[i-1][0] = (i-1 == 0) ? 0 : mu;
		B[i-1][1] = 2;
		B[i-1][2] = (i-1 == n-2) ? 0 : lambda;
	}

	dVector coefs = solTriD_Thomas(B,d,sqrt(getMachineEpsilon<double>()));
	
	std::vector<double> M;
	M.push_back(0.0);

	for (int i = 0; i < n-1;i++)
		M.push_back(coefs[i][0]);

	M.push_back(0.0);
	return M;
}

// - - - - - - - - Norms - - - - - - - - 

double normL1(const dVector& x)
{
	// LOG("normL1()");
	double sum = 0;
	int rows, cols;
	rows = x.getRows();
	cols = x.getCols();

	if (rows == 1)
		for (int i = 0; i < cols; i++) sum += std::abs(x[0][i]);
	else if (cols == 1)
		for (int i = 0; i < rows; i++) sum += std::abs(x[i][0]);

	return sum;
}

double normL2(const dVector& x)
{
	// LOG("normL2()");
	double sum = 0;
	int rows, cols;
	rows = x.getRows();
	cols = x.getCols();

	if (rows == 1)
		for (int i = 0; i < cols; i++) sum += pow(x[0][i],2);
	else if (cols == 1)
		for (int i = 0; i < rows; i++) sum += pow(x[i][0],2);

	return sqrt(sum);
}

double normFrobenius(const dMatrix& A)
{
	double n = A.getRows(), m = A.getCols();

	double sumout = 0;
	for (int i = 0; i < n; i++)
	{
		double sumin = 0;
		for (int j = 0; j < m; j++)
		{
			sumin += std::abs(A[i][j])*std::abs(A[i][j]);
		}
		sumout += sumin;
	}

	return sqrt(sumout);
	
}

// - - - - - - - - Error calculations - - - - - - - - 

double MSE(const dVector& Yreal, const dVector& Yhat)
{
	if (Yreal.getRows() != Yhat.getRows() || Yreal.getCols() != Yhat.getCols())
	{
		perror("Vector size does not match");
		return std::numeric_limits<double>::infinity();
	}

	double error = 0;

	if (Yreal.getRows() == 1)
	{
		int n = Yreal.getCols();
		for (int i = 0; i<n; i++)
			error += pow(Yhat[0][i]-Yreal[0][i],2);

		error/=n;
	}

	else if (Yreal.getCols() == 1)
	{
		int n = Yreal.getRows();
		for (int i = 0; i<n; i++)
			error += pow(Yhat[i][0]-Yreal[i][0],2);

		error/=n;
	}

	return error;
}
