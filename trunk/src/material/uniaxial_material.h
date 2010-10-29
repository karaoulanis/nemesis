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

#ifndef _UNIAXIALMATERIAL_H
#define _UNIAXIALMATERIAL_H

#include "material/material.h"

/**
 * The Uniaxial Material Class.                                                
 */
class UniaxialMaterial: public Material
{
protected:
	double sTrial;
	double sConvg;
	double eTotal;
public:
	UniaxialMaterial();
	UniaxialMaterial(int ID,double rho,double aT);
	~UniaxialMaterial();

	virtual UniaxialMaterial* getClone()=0;
	virtual void setStrain(const double De)=0;
	virtual double getC()=0;
	inline void setStress(const double s)	{sTrial=s;		}
	inline void addStress(const double s)	{sTrial+=s;		}///@todo: check
	inline double getStress()			{return sTrial;	}
	
	// Tracker member functions
	void track();
};

#endif
