
#include "material/hardening.h"

Hardening::Hardening()
{
}
double Hardening::geth(const Vector& v)
{
	return sqrt(2./3.*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));
}
const Vector& Hardening::gethds(const Vector& /*sigma*/,const double /*kappa*/)
{
	static Vector a(3,0.);
	return a;
}
double Hardening::gethdk(const Vector& /*sigma*/,const double /*kappa*/)
{
	return 0;
}

