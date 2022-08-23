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

#ifndef functional
	#include <functional>
#endif

#ifndef fstream
	#include <fstream>
	#include <sstream>
#endif

//integral from [a,b] of f ---> (f,[a,b])
double RTrapezoidal(std::function<double(double)>f, double a, double b, int n); //R_{f;[a,b]}(n,0) <--- this is a terrible idea
dMatrix getRombergMat(std::function<double(double)> f, double a, double b, int n);
