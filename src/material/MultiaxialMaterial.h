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

#ifndef _MULTIAXIALMATERIAL_H
#define _MULTIAXIALMATERIAL_H

#include "material/Material.h"
#include "numeric/matrix.h"
#include "numeric/vector.h"

/**
 * The Multiaxial Material Class.                                                
 */
class MultiaxialMaterial: public Material
{
protected:
	Vector sTrial;
	Vector sConvg;
	Vector eTrial;
	Vector eTotal;
	static Matrix C;
public:
	MultiaxialMaterial();
	MultiaxialMaterial(int ID,double rho,double aT);
	~MultiaxialMaterial();

	virtual MultiaxialMaterial* getClone()=0;
	virtual void setStrain(const Vector& De)=0;
	virtual const Matrix& getC()=0;
	void  setStress(const Vector& s){sTrial =s; sConvg =s;}
	void  addStress(const Vector& s){sTrial+=s; sConvg+=s;}
	const Vector& getStress()		{return sTrial;}
	virtual bool isPlastic()		{return false;}
};
#endif
