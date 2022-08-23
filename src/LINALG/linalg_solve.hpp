#ifndef MAT_TRUNCATION
	#include "../MATRIX/matrix.hpp"
#endif

#ifndef LEPSILON
	#include "../MATRIX/epsilon.hpp"
#endif

#ifndef cmath
	#include <cmath>
#endif

#ifndef vector
	#include <vector>
#endif

#ifndef algorithm
	#include <algorithm>
#endif

#ifndef limits
	#include <limits>
#endif

#ifndef fstream
	#include <fstream>
	#include <sstream>
#endif

//A = arbitratry matrix
//L = lower triangular matrix
//U = upper triangular matrix
//b = solution vector
//p = pertutation vector

//tau = tolerance

//Substitutions
dVector forwardSubstitution(const dMatrix&, const dVector&, double); //A, b, tau
dVector backwardSubstitution(const dMatrix&, const dVector&, double); //A, b, tau

//Solvers
//-explicit
dVector genSolLU (const dMatrix&, const dMatrix&, const dVector&, unsigned int*, double); //L, U, b, p, tau
dVector solFactLU (const dMatrix, const dVector); //A, b 
dVector solFactChol (const dMatrix, const dVector, double);
dVector solTriD_Thomas (const dMatrix&, const dVector&, double);
dVector fitLeastSquares (const dMatrix&, const dVector&);
//-iterative
int solTriD_GS (const dMatrix&, const dVector&, dVector&, unsigned int, double, double&);
dVector approxSolSVD(const dVector&, const dMatrix&, const dMatrix&, std::vector<double>&, unsigned int);

//Factorizations
bool factLU(const dMatrix, dMatrix&, dMatrix&, unsigned int*, double); //A, L, U, p, tau
dMatrix factChol(const dMatrix, double);
std::vector<double> iterativeSVD (dMatrix, double, unsigned int, dMatrix&, dMatrix&);

//Helper functions
dMatrix reduceTriDMatrix(const dMatrix&);
dMatrix getDataSet (const char[]);
dMatrix generatePhi (int degree, const dMatrix& data);
double evaluateHorner(int, double*, double);
double evaluateSpline (double x, const dMatrix& data, const std::vector<double> M);
void exportPolinomial(int, double*, int, double*, const std::string);
void exportSpline (double a, double b, const dMatrix& data, const std::vector<double> M, const std::string& outPath);
double matrixConditionNumber(const dMatrix&);
std::vector<double> getD2CSplineCoef(const dMatrix&);

//Norms
double normL1(const dVector&);
double normL2(const dVector&);
double normFrobenius(const dMatrix&);

//Error calculations
double MSE(const dVector&, const dVector&);
