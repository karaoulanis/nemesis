/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2007 F.E.Karaoulanis [http://www.nemesis-project.org]  *
*                                                                             *
*   This program is free software; you can redistribute it and/or modify      *
*   it under the terms of the GNU General Public License version 3, as        *
*   published by the Free Software Foundation.                                *
*                                                                             *
*   This program is distributed in the hope that it will be useful,           *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*   GNU General Public License for more details.                              *
*                                                                             *
*   You should have received a copy of the GNU General Public License         *
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
******************************************************************************/

//*****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
//*****************************************************************************

#ifndef _NUMERIC_H
#define _NUMERIC_H

#include <cmath>
#include <limits>
#include <cstdio>
using namespace std;

namespace num
{
	const double pi=4.*atan(1.0);
	const double sq2=sqrt(2.);
	const double sq3=sqrt(3.);
	const double sq6=sqrt(6.);
	const double d13=1./3.;
	const double d23=2./3.;
	const double d43=4./3.;
	const double d16=1./6.;
	const double eps=numeric_limits<double>::epsilon();
	inline bool equal(double d1,double d2)
	{
		return fabs(d2-d1)<eps	? true : false;
	}
	inline bool smaller(double d1,double e)
	{
		return fabs(d1)<e		? true : false;
	}
	inline bool smaller(double d1,double d2,double e)
	{
		return fabs(d1-d2)<e	? true : false;
	}
	inline bool tiny(double d)
	{
		return fabs(d)<eps		? true : false;
	}
	inline double deg2rad(double d)
	{
		return d*pi/180.;
	}
	inline double rad2deg(double d)
	{
		return d*180./pi;
	}
	inline void print_d(double d,int total,int decimal)
	{
		char format[64];
		sprintf(format,"%% %d.%df",total,decimal);
		int fw=total-1-1-decimal;
		if(fabs(d)<pow(10.0,(double)fw))	printf(format,d);
		else for(int i=0;i<total;i++)		printf("*");
	}
	inline void print_i(int n,int total)
	{
		char format[64];
		sprintf(format,"%% %dd",total);
		int fw=total-1;
		///@todo Check for abs(int)
		if(fabs((double)n)<pow(10.0,(double)fw))	printf(format,n);
		else for(int i=0;i<total;i++)		printf("*");
	}
	inline double sign(double d)
	{
		return d<0	? -1 : 1;
	}
}
#endif
