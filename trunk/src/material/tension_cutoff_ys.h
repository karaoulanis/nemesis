
#ifndef _TENSIONCUTOFFYS_H
#define _TENSIONCUTOFFYS_H

#include "material/ys.h"

class TensionCutOffYS: public YS
{
private:
	double T;
public:
	TensionCutOffYS(double T_);
	double getf(const Vector& sigma,const double kappa);
	const Vector& getdfds(const Vector& sigma,const double kappa);
	const Matrix& getd2fdsds(const Vector& sigma,const double kappa);
	double getdfdk(const Vector& sigma,const double kappa);
	const Vector& getf2dkds(const Vector& sigma,const double kappa);
};

#endif
