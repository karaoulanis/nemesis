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

#include <VM.h>

VM::VM()
{
}
VM::VM(double s0_,double K_)
{
	s0=s0_;
	K=K_;
}
VM::~VM()
{
}
double VM::get_f(const Vector& s,const Vector& e)
{
	Vector r=this->get_dfds(s,e);
	double eq=sqrt(2/3.)*e.twonorm();
	//cout<<num::sq3*sqrt(s.J2())-(s0+K*eq)<<endl;
	return num::sq3*sqrt(s.J2())-(s0+K*eq);
}
void VM::find_C(const Vector& s,const Vector& e)
{
	C1=0.;
	C2=sqrt(0.75)/sqrt(s.J2());
	C3=0.;
	C11=0.;
	C22=-0.25*num::sq3*pow(s.J2(),-1.5);
	C23=0.; C32=0.; C33=0.;
}
const double  VM::get_dfdq(const Vector& s,const Vector& e)
{
	return -1;
}
const Vector& VM::get_df2dsq(const Vector& s,const Vector& e)
{
	static Vector ret(6,0.);
	ret.clear();
    return ret;
}
const double  VM::get_df2dqq(const Vector& s,const Vector& e)
{
	return 0.;
}
const double  VM::get_H(const Vector& s,const Vector& e)
{
	Vector r=this->get_dfds(s,e);
	return K;
//	return -K*sqrt(2/3.)*sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]+2*r[3]*r[3]+2*r[4]*r[4]+2*r[5]*r[5]);
}
