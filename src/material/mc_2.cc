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

#include "material/MC_2.h"

MC_2::MC_2()
{
}
MC_2::MC_2(double c_,double phi_)
:MC(c_,phi_)
{
}
MC_2::~MC_2()
{
}
void MC_2::find_A(const Vector& s,double& A,double& dA,double& d2A)
{
	double theta=s.theta();
	A  = 0.5*cos(theta)*(1.+sin(phi))+0.5*sin(theta)/num::sq3*(sin(phi)-3.);
	dA =-0.5*sin(theta)*(1.+sin(phi))+0.5*cos(theta)/num::sq3*(sin(phi)-3.);
	d2A=-0.5*cos(theta)*(1.+sin(phi))-0.5*sin(theta)/num::sq3*(sin(phi)-3.);
}
