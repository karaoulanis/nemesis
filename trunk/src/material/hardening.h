
#ifndef _HARDENING_H
#define _HARDENING_H

#include "numeric/matrix.h"

class Hardening
{
private:
public:
	Hardening();
	double geth(const Vector& v);
	const Vector& gethds(const Vector& sigma,const double kappa);
	const double gethdk(const Vector& sigma,const double kappa);
};

#endif