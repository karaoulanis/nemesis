/*******************************************************************************
* nemesis. an experimental finite element code.                                *
* Copyright (C) 2004-2010 F.E.Karaoulanis [http://www.nemesis-project.org]     *
*                                                                              *
* This program is free software; you can redistribute it and/or modify         *
* it under the terms of the GNU General Public License version 3, as           *
* published by the Free Software Foundation.                                   *
*                                                                              *
* This program is distributed in the hope that it will be useful,              *
* but WITHOUT ANY WARRANTY; without even the implied warranty of               *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
* GNU General Public License for more details.                                 *
*                                                                              *
* You should have received a copy of the GNU General Public License            *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
// *****************************************************************************

#include "material/mc_0.h"

MC_0::MC_0()
{
}
MC_0::MC_0(double c_,double phi_)
:MC(c_,phi_)
{
}
MC_0::~MC_0()
{
}
void MC_0::find_A(const Vector& s,double& A,double& dA,double& d2A)
{
	double theta=s.theta();
	A  = cos(theta)-sin(theta)*sin(phi)/num::sq3;
	dA =-sin(theta)-cos(theta)*sin(phi)/num::sq3;
	d2A=-cos(theta)+sin(theta)*sin(phi)/num::sq3;
}
