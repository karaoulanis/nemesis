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
VM::VM(double s0_)
{
	s0=s0_;
}
VM::~VM()
{
}
double VM::get_f(const Vector& s,const Vector& e)
{
	return num::sq3*sqrt(s.J2())-s0;
}
void VM::find_C(const Vector& s,const Vector& e)
{
	C1=0.;
	C2=sqrt(0.75)/sqrt(s.J2());
	C3=0.;
	C11=0.;
	C22=-0.25*num::sq3/pow(s.J2(),-1.5);
	C23=0.; C32=0.; C33=0.;
}
