
#ifndef _YS_H
#define _YS_H

#include <vector>

#include "main/nemesis_debug.h"
#include "numeric/matrix.h"

class YS
{
protected:
	double s1,s2,s3;
	double I1,J2;
	Vector a,a1,a2;
	Matrix da,da2,da11,da22;
	bool active;
public:
	YS();
	void setSigma(const Vector& s);

	// functions to be overwritten
	virtual double getf(const Vector& sigma,const double kappa)=0;
	virtual const Vector& getdfds(const Vector& sigma,const double kappa)=0;
	virtual const Matrix& getd2fdsds(const Vector& sigma,const double kappa)=0;
	virtual double getdfdk(const Vector& sigma,const double kappa)=0;
	virtual const Vector& getf2dkds(const Vector& sigma,const double kappa)=0;
};
#endif
