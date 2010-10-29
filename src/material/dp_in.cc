/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2010 F.E.Karaoulanis [http://www.nemesis-project.org]  *
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

#include "material/dp_in.h"

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
	double Kf=0.2;
	//cout<<phi<<'\t'<<Kf*a<<endl;
	double Kc= 0.;
	double D= 3*2*sin(phi+Kf*a)/(num::sq3*(3+sin(phi+Kf*a)));
	double so=6*(c+Kc*a)*cos(phi-Kf*a)/(num::sq3*(3+sin(phi-Kf*a)));
	return D*s.I1()+sqrt(s.J2())-so;
}
void DP_in::find_C(const Vector& s,const double a)
{
	double Kf=0.;
	C1=2*sin(phi+Kf*a)/(num::sq3*(3+sin(phi+Kf*a)));
	//C1=2*sin(phi)/(num::sq3*(3+sin(phi)));
	C2=0.5/sqrt(s.J2());
	C3=0.;
	C11=0.;
	C22=-0.25*pow(s.J2(),-1.5);
	C23=0.; C32=0.; C33=0.;
}
double DP_in::get_dfda(const Vector& s,const double a)
{
	double Kf=0.;
/*
	double Kf=0.0075;
	double sf=3+sin(phi+Kf*a);
	double cf=cos(phi+Kf*a);
	double d1=-6/(sqrt(3.)*sf*sf)*(cf*s.I1()+c*(3*cf+1))*Kf;
	double d2=+6*cf/(sqrt(3.)*sf)*Kc;

	return -(d1+d2);
	//return 0;
*/
	return -2*Kf*sqrt(3.)*(cos(phi+Kf*a)*s.I1()+3*c*sin(phi+Kf*a)+c)/
		(-10-6.*sin(phi+Kf*a)+pow(cos(phi+Kf*a),2));
}
const Vector& DP_in::get_df2dsa(const Vector& /*s*/,const double a)
{
	static Vector ret(6,0.);
	ret.clear();
	double Kf=0.;
	double d=-2*Kf*sqrt(3.)*cos(phi+Kf*a)/(-10.-6.*sin(phi+Kf*a)+pow(cos(phi+Kf*a),2));
	ret[0]=d;
	ret[1]=d;
	ret[2]=d;

/*
	double Kc=0.;
	double Kf=0.01;
	double sf=3.+sin(phi);
	double cf=cos(phi);
	ret[0]=(6./(sqrt(3.)*sf*sf))*cf*Kf;
	ret[1]=(6./(sqrt(3.)*sf*sf))*cf*Kf;
	ret[2]=(6./(sqrt(3.)*sf*sf))*cf*Kf;
*/
    return ret;
}
double DP_in::get_df2daa(const Vector& /*s*/,const double /*a*/)
{
	return 0.;
}
