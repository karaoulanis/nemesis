/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2007 F.E.Karaoulanis [http://www.nemesis-project.org]  *
*                                                                             *
*   This program is free software; you can redistribute it and/or modify      *
*   it under the terms of the GNU General Public License version 2, as        *
*   published by the Free Software Foundation.                                *
*                                                                             *
*   This program is distributed in the hope that it will be useful,           *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*   GNU General Public License for more details.                              *
*                                                                             *
*   You should have received a copy of the GNU General Public License along   *
*   with this program; if not, write to the Free Software Foundation, Inc.,   *
*   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.               *
******************************************************************************/

//*****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
//*****************************************************************************

#include <DP_in.h>

DP_in::DP_in()
{
}
DP_in::DP_in(double c_,double phi_)
{
	c=c_;
	phi=phi_;
}
DP_in::~DP_in()
{
}
double DP_in::get_f(const Vector& s,const double a)
{
	double Kf=10.;
	double D=   2*sin(phi+Kf*a)/(num::sq3*(3+sin(phi+Kf*a)));
	double so=6*c*cos(phi+Kf*a)/(num::sq3*(3+sin(phi+Kf*a)));
	return D*s.I1()+sqrt(s.J2())-so;
}
void DP_in::find_C(const Vector& s,const double a)
{
	double Kf=10.;
	C1=2*sin(phi+Kf*a)/(num::sq3*(3+sin(phi+Kf*a)));
	//C1=2*sin(phi)/(num::sq3*(3+sin(phi)));
	C2=0.5/sqrt(s.J2());
	C3=0.;
	C11=0.;
	C22=-0.25*pow(s.J2(),-1.5);
	C23=0.; C32=0.; C33=0.;
}
const double  DP_in::get_dfda(const Vector& s,const double a)
{
	double Kc=10.;
	double Kf=0.0075;
	double sf=3+sin(phi+Kf*a);
	double cf=cos(phi+Kf*a);
	double d1=-6/(sqrt(3.)*sf*sf)*(cf*s.I1()+c*(3*cf+1))*Kf;
	double d2=+6*cf/(sqrt(3.)*sf)*Kc;

	return -(d1+d2);
	//return 0;
}
const Vector& DP_in::get_df2dsa(const Vector& s,const double a)
{
	static Vector ret(6,0.);
	ret.clear();

	double Kc=0.;
	double Kf=10.;
	double sf=3.+sin(phi);
	double cf=cos(phi);
	ret[0]=(6./(sqrt(3.)*sf*sf))*cf*Kf;
	ret[1]=(6./(sqrt(3.)*sf*sf))*cf*Kf;
	ret[2]=(6./(sqrt(3.)*sf*sf))*cf*Kf;

    return ret;
}
const double  DP_in::get_df2daa(const Vector& s,const double a)
{
	return 0.;
}
