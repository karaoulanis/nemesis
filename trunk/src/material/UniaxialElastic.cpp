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

#include <UniaxialElastic.h>

UniaxialElastic::UniaxialElastic()
{
}
UniaxialElastic::UniaxialElastic(int ID,double E,double nu,double rho,double aT)
:UniaxialMaterial(ID,rho,aT)
{
	// Material parameters
	MatParams[0]=E;
	MatParams[1]=nu;
	myTag=TAG_MATERIAL_UNIAXIAL_ELASTIC;
}
UniaxialMaterial* UniaxialElastic::getClone()
{
	// Material parameters
	double E   =MatParams[ 0];
	double nu  =MatParams[ 1];
	double rho =MatParams[30];
	double aT  =MatParams[31];
	// Create clone and return
	UniaxialMaterial* clone=new UniaxialElastic(myID,E,nu,rho,aT);
	return clone;
}
void UniaxialElastic::setStrain(const double de)
{
	sTrial=sConvg+MatParams[0]*de;
}
const double UniaxialElastic::getC()
{
	return MatParams[0];
}
void UniaxialElastic::commit()
{
	sConvg=sTrial;
}
