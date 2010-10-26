
#ifndef _DRUCKERPRAGERYS_H
#define _DRUCKERPRAGERYS_H

#include "material/ys.h"

class DruckerPragerYS: public YS
{
private:
	double c0,phi0,Kc,Kphi;
public:
	DruckerPragerYS(double c_,double phi_,double Kc_,double Kphi_);
	const double getf(const Vector& sigma,const double kappa);
	const Vector& getdfds(const Vector& sigma,const double kappa);
	const Matrix& getd2fdsds(const Vector& sigma,const double kappa);
	const double getdfdk(const Vector& sigma,const double kappa);
	const Vector& getf2dkds(const Vector& sigma,const double kappa);
};

#endif
