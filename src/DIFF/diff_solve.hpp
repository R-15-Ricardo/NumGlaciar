#pragma once

#include "../MATRIX/matrix.hpp"
#include "../MATRIX/epsilon.hpp"
#include "../LINALG/linalg_solve.hpp"
#include <functional>

#include <string>
#include <fstream>

bool nanVector (dVector x);

double derivative(std::function<double(double)> f, double x, double h);
dVector grad(std::function<double(dVector)> f, const dVector& x, double h);

dVector autograd(std::function<double(dVector)> f, dVector start, double stepSize, double tolerance, int verbose);
dVector optimize(std::function<dVector(dVector)> df, dVector start,  double stepSize, double tolerance, int verbose);

void plot3d(const std::string outPath, std::function<double(dVector)> f, double xmin, double xmax, double ymin, double ymax, int finess = 100);

