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

#include "material/tc.h"

TC::TC()
{
}
TC::TC(double t_)
{
	t=t_;
}
TC::~TC()
{
}
double TC::get_f(const Vector& s,const double q)
{
	return t+s.I1();
}
const Vector& TC::get_dfds(const Vector& s,const double a)
{
	myVector[0]=1.0;
	myVector[1]=1.0;
	myVector[2]=1.0;
	myVector[3]=0.0;
	myVector[4]=0.0;
	myVector[5]=0.0;
	return myVector;
}
