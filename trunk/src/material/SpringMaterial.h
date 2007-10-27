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

#ifndef _SPRINGMATERIAL_H
#define _SPRINGMATERIAL_H

#include <Material.h>

/**
 * The Single Dof Material Class.                                                
 */
class SpringMaterial: public Material
{
protected:
	Vector sTrial;
	Vector sConvg;
	Vector eTrial;
	Vector eTotal;
	Vector Ct;
	int nDim;
public:
	SpringMaterial();
	SpringMaterial(int ID);

	// Get clone
	virtual SpringMaterial* getClone()=0;
	virtual void setStrain(const Vector& De)=0;

	const Vector& getC();
	void commit();
	inline void setStress(const Vector& s)	{sTrial=s;		}
	inline void addStress(const Vector& s)	{sTrial+=s;		}
	inline const Vector& getStress()		{return sTrial;	}

	// Tracker member functions
	void track();
};

#endif
