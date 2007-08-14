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

#include <Crack.h>

Crack::Crack()
{
}
Crack::Crack(int ID,double xS,double yS,double xT,double yT)
{
	myCrackPoints.push_back(new CrackPoint(xS,yS));
	myCrackPoints.push_back(new CrackPoint(xT,yT));
}
Crack::~Crack()
{
	Containers::vector_delete(myCrackPoints);
}
const vector<CrackPoint*>& Crack::getCrack()
{
	return myCrackPoints;
}
CrackPoint* Crack::getCrackPoint(int n)
{
	return myCrackPoints[n];
}