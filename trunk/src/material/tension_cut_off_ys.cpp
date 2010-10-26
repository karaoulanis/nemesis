
#include "material/TensionCutOffYS.h"


TensionCutOffYS::TensionCutOffYS(double T_)
{
	T=T_;
}
const double TensionCutOffYS::getf(const Vector& sigma,const double kappa)
{
	this->setSigma(sigma);
	return I1-T;
}
const Vector& TensionCutOffYS::getdfds(const Vector& sigma,const double kappa)
{
	this->setSigma(sigma);
	a=a1;
	return a;
}
const Matrix& TensionCutOffYS::getd2fdsds(const Vector& sigma,const double kappa)
{
	da.clear();
	return da;
}
const double TensionCutOffYS::getdfdk(const Vector& sigma,const double kappa)
{
	return 0;
}
const Vector& TensionCutOffYS::getf2dkds(const Vector& sigma,const double kappa)
{
	a.clear();
	return a;
}

