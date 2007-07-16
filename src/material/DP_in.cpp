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
DP_in::DP_in(double s0_)
{
	s0=s0_;
}
DP_in::~DP_in()
{
}
double DP_in::get_f(const Vector& s)
{
	return num::sq3*sqrt(s.J2())-s0;
}
void DP_in::find_C(const Vector& s)
{
	C1=0.;
	C2=sqrt(0.75)/sqrt(s.J2());
	C3=0.;
	C22=-0.25*num::sq3/pow(s.J2(),-1.5);
	C23=0.; C32=0.; C33=0.;
}
