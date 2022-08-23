#include "diff_solve.hpp"

bool nanVector (dVector x)
{
	int n = x.getRows();
	for (int i = 0; i < n; i++)
	{
		if (std::isnan(x[i][0]))
			return true;
	}
	return false;
}

double derivative(std::function<double(double)> f, double x, double h)
{
	double back = f(x+h);
	double forth = f(x-h);

	return (back - forth)/(2*h);
}

dVector grad(std::function<double(dVector)> f, const dVector& x, double h)
{
	int n = x.getRows();

	dVector gradVector (n,1,dVector::MatrixType::zeros);
	for (int i = 0; i < n; i++)
	{
		auto f_i = [=](double x_0)->double {
			dVector trunc =  x;
			trunc[i][0] = x_0;
			return f(trunc);
		};
		gradVector[i][0] = derivative(f_i,x[i][0],h);
	}

	return gradVector;
}

dVector autograd(std::function<double(dVector)> f, dVector start, double stepSize, double tolerance, int verbose)
{
	dVector offset(start.getRows(),1,dVector::MatrixType::ones);
	dVector stepNow = 10*offset + start;
	dVector stepLast = start;

	std::fstream file;
	if (verbose == 2)
		file.open("gradient_path.txt",std::ios::out);


	double diff = normL2(stepNow - stepLast);

	while (diff > tolerance)
	{
		stepNow = stepLast - stepSize*grad(f,stepLast,tolerance);
		diff = normL2(stepNow - stepLast);
		if (nanVector(stepNow))
			break;

		if (verbose == 1)
			std::cout<<"argmin now: "<<stepNow.T()<<std::endl;

		if (verbose == 2)
		{
			double p1x = stepLast[0][0], p1y = stepLast[1][0], p1z = f(stepLast);
			double p2x = stepNow[0][0], p2y = stepNow[1][0], p2z = f(stepNow);

			file<<p1x<<" "<<p1y<<" "<<p1z<<" ";
			file<<p2x - p1x<<" "<<p2y - p1y<<" "<<p2z - p1z<<std::endl;
		}

		stepLast = stepNow;

	}

	if (verbose == 2)
		file.close();

	return stepNow;
}

dVector optimize(std::function<dVector (dVector)> df, dVector start,  double stepSize, double tolerance, int verbose)
{
	dVector offset(start.getRows(),1,dVector::MatrixType::ones);
	dVector stepNow = 10*offset + start;
	dVector stepLast = start;

	double diff = normL2(stepNow - stepLast);


	while (diff > tolerance)
	{
		stepNow = stepLast - stepSize*df(stepLast);
		diff = normL2(stepNow - stepLast);
		if (nanVector(stepNow))
			break;

		if (verbose == 1)
			std::cout<<"argmin now: "<<stepNow.T()<<std::endl;

		if (verbose == 2)
		{

		}

		stepLast = stepNow;
	}

	return stepNow;
}

void plot3d(const std::string outPath, std::function<double(dVector)> f, double xmin, double xmax, double ymin, double ymax, int finess)
{
	using namespace std;
	fstream file;
	file.open(outPath.c_str(),ios::out);

	if (file.is_open())
	{
		double stepx = (xmax-xmin)/finess;
		double stepy = (ymax-ymin)/finess;

		for (int xi = 0; xi < finess; xi++)
		{
			for(int yi = 0; yi < finess; yi++)
			{
				double xval = xmin + xi*stepx, yval = ymin + yi*stepy;
				dVector p(2,1);

				p[0][0] = xval;
				p[1][0] = yval;

				file<<xval<<" "<<yval<<" "<<f(p)<<std::endl;
			}
			file<<std::endl;
		}
	}
	file.close();
}